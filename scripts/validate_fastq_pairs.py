#!/usr/bin/env python3
"""Standalone paired-FASTQ validator for fastBact Phase F26 (F26-DATA-002).

Validates Illumina paired-end FASTQ gzip archives produced by the F26
download pipeline (scripts/download_f26_datasets.sh):

  a) file existence and readability
  b) non-empty gzip integrity (full stream decompression; CRC/truncation)
  c) exact line count divisibility by 4
  d) exact spot/pair count equality between R1 and R2
  e) basic read ID pair synchronization between mates
     (first header token compared after stripping optional /1 and /2
     suffixes, so both the ``@id/1``-``@id/2`` and the fasterq-dump
     ``@id``-``@id`` conventions are accepted)
  f) structured JSON report of all findings

Optionally cross-checks measured R1 read counts against the frozen
``F26_DATASET_MANIFEST.tsv`` ``total_spots`` column (--manifest). By
default a manifest mismatch is reported but does not fail validation;
pass --enforce-manifest-spots to make it a failure.

Exit codes: 0 = all accessions PASS, 1 = at least one FAIL, 2 = usage or
environment error. Memory usage is O(1): files are streamed, never loaded.

Usage:
    python3 scripts/validate_fastq_pairs.py \
        --fastq-dir /path/to/fastq \
        --manifest F26_DATASET_MANIFEST.tsv \
        --output f26_validation_report.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import zlib
from datetime import datetime, timezone
from typing import Iterator, Optional

SCHEMA_VERSION = "1.0"
TOOL_NAME = "validate_fastq_pairs.py"

# Header positions inside a 4-line FASTQ record (0-based).
HEADER_OFFSET = 0


class FastqStreamError(Exception):
    """Raised when a gzip FASTQ stream cannot be fully decompressed."""


def _open_text_stream(path: str) -> io.TextIOBase:
    """Open a (possibly gzip-compressed) FASTQ for streaming text reads.

    Prefers an external ``pigz -dc`` pipe for speed; falls back to the
    stdlib gzip module when pigz is unavailable or explicitly disabled
    via FASTQ_VALIDATOR_NO_PIGZ=1.
    """
    if os.environ.get("FASTQ_VALIDATOR_NO_PIGZ") != "1" and shutil.which("pigz"):
        proc = subprocess.Popen(
            ["pigz", "-dc", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        wrapper = io.TextIOWrapper(proc.stdout, encoding="utf-8", errors="replace")

        class PigzTextStream:
            """Text stream that fails closed if the pigz process failed."""

            def __init__(self) -> None:
                self._proc = proc
                self._closed = False

            def __iter__(self) -> Iterator[str]:
                return self

            def __next__(self) -> str:
                try:
                    line = next(wrapper)
                except StopIteration:
                    self._finish()
                    raise
                return line

            def _finish(self) -> None:
                if self._closed:
                    return
                self._closed = True
                wrapper.close()
                rc = self._proc.wait()
                if rc != 0:
                    raise FastqStreamError(
                        f"pigz -dc exited with status {rc} (corrupt or "
                        f"truncated gzip stream): {path}"
                    )

            def close(self) -> None:
                if not self._closed:
                    try:
                        # Drain remaining output so pigz is not killed by
                        # SIGPIPE before its exit status can be observed.
                        for _ in self:
                            pass
                    except (FastqStreamError, OSError):
                        pass
                self._proc.stdout.close()
                if self._proc.poll() is None:
                    self._proc.terminate()
                    self._proc.wait()

        return PigzTextStream()  # type: ignore[return-value]

    try:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    except OSError as exc:  # unreadable / not gzip at all
        raise FastqStreamError(f"cannot open gzip stream: {path}: {exc}") from exc


def _normalized_read_id(header_line: str) -> str:
    """Extract the read identifier from a FASTQ header line.

    Takes the first whitespace-delimited token, drops the leading ``@``
    and an optional trailing ``/1`` or ``/2`` mate suffix.
    """
    token = header_line.rstrip("\r\n").split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def inspect_file(path: str) -> dict:
    """Stream one FASTQ file and collect structural integrity metrics."""
    info: dict = {
        "path": path,
        "exists": os.path.isfile(path),
        "readable": False,
        "size_bytes": None,
        "gzip_integrity_ok": False,
        "line_count": None,
        "lines_divisible_by_4": None,
        "read_count": None,
        "malformed_records": 0,
        "error": None,
    }
    if not info["exists"]:
        info["error"] = "file does not exist"
        return info
    info["size_bytes"] = os.path.getsize(path)
    if not os.access(path, os.R_OK):
        info["error"] = "file is not readable (permission)"
        return info
    if info["size_bytes"] == 0:
        info["error"] = "file is empty (0 bytes)"
        return info
    info["readable"] = True

    line_count = 0
    try:
        stream = _open_text_stream(path)
        try:
            for line in stream:
                if line_count % 4 == HEADER_OFFSET and not line.startswith("@"):
                    info["malformed_records"] += 1
                line_count += 1
        finally:
            stream.close()
    except (FastqStreamError, OSError, EOFError, zlib.error) as exc:
        info["error"] = f"gzip/decompression failure after {line_count} lines: {exc}"
        return info

    info["gzip_integrity_ok"] = True
    info["line_count"] = line_count
    info["lines_divisible_by_4"] = line_count % 4 == 0
    if info["lines_divisible_by_4"]:
        info["read_count"] = line_count // 4
    else:
        info["error"] = f"line count {line_count} is not divisible by 4"
    return info


def iter_read_ids(path: str) -> Iterator[str]:
    """Yield the normalized read ID of every record header, streaming."""
    stream = _open_text_stream(path)
    try:
        for line_no, line in enumerate(stream):
            if line_no % 4 == HEADER_OFFSET:
                yield _normalized_read_id(line)
    finally:
        stream.close()


def check_id_sync(r1_path: str, r2_path: str, max_examples: int) -> dict:
    """Compare mate read IDs pairwise; O(1) memory via concurrent streams."""
    result: dict = {
        "checked": True,
        "mismatch_count": 0,
        "compared_pairs": 0,
        "examples": [],
        "error": None,
    }
    try:
        ids1 = iter_read_ids(r1_path)
        ids2 = iter_read_ids(r2_path)
        for idx, (id1, id2) in enumerate(zip(ids1, ids2)):
            result["compared_pairs"] = idx + 1
            if id1 != id2:
                result["mismatch_count"] += 1
                if len(result["examples"]) < max_examples:
                    result["examples"].append(
                        {"pair_index": idx, "r1_id": id1, "r2_id": id2}
                    )
    except (FastqStreamError, OSError, EOFError, zlib.error) as exc:
        result["error"] = str(exc)
    return result


def load_manifest(path: str) -> dict:
    """Map SRA accession -> {dataset_id, total_spots} from the F26 manifest."""
    table: dict = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("sra_accession") or "").strip()
            if accession:
                table[accession] = {
                    "dataset_id": (row.get("dataset_id") or "").strip(),
                    "total_spots": (row.get("total_spots") or "").strip(),
                }
    return table


def find_accessions(fastq_dir: str) -> list[str]:
    """Auto-discover ``<acc>_1.fastq.gz`` accessions in a directory."""
    found = set()
    for name in sorted(os.listdir(fastq_dir)):
        if name.endswith("_1.fastq.gz"):
            found.add(name[: -len("_1.fastq.gz")])
    return sorted(found)


def validate_accession(accession: str, fastq_dir: str, manifest: Optional[dict],
                       max_examples: int, enforce_manifest: bool) -> dict:
    """Run every check for one accession and return its report entry."""
    r1_path = os.path.join(fastq_dir, f"{accession}_1.fastq.gz")
    r2_path = os.path.join(fastq_dir, f"{accession}_2.fastq.gz")
    findings: list[str] = []

    r1 = inspect_file(r1_path)
    r2 = inspect_file(r2_path)

    for mate, info in (("R1", r1), ("R2", r2)):
        if info["error"]:
            findings.append(f"{mate}: {info['error']}")

    read_count_equal: Optional[bool] = None
    if r1["read_count"] is not None and r2["read_count"] is not None:
        read_count_equal = r1["read_count"] == r2["read_count"]
        if not read_count_equal:
            findings.append(
                f"pair count mismatch: R1={r1['read_count']} != "
                f"R2={r2['read_count']}"
            )

    id_sync = {
        "checked": False,
        "mismatch_count": None,
        "compared_pairs": None,
        "examples": [],
        "error": None,
    }
    if (
        r1["gzip_integrity_ok"]
        and r2["gzip_integrity_ok"]
        and read_count_equal
    ):
        id_sync = check_id_sync(r1_path, r2_path, max_examples)
        id_sync["checked"] = True
        if id_sync["error"]:
            findings.append(f"read ID synchronization: {id_sync['error']}")
        elif id_sync["mismatch_count"]:
            findings.append(
                f"read ID desynchronization in {id_sync['mismatch_count']} "
                f"of {id_sync['compared_pairs']} pairs"
            )
    else:
        findings.append("read ID synchronization skipped (files not both valid)")

    # Singleton file produced by fasterq-dump --split-3 (informational).
    singleton = None
    for cand in (f"{accession}.fastq.gz", f"{accession}.fastq"):
        path = os.path.join(fastq_dir, cand)
        if os.path.isfile(path):
            singleton = {"path": path, "size_bytes": os.path.getsize(path)}
            break

    manifest_info = None
    if manifest is not None:
        entry = manifest.get(accession)
        if entry is None:
            manifest_info = {"present_in_manifest": False}
            findings.append("accession not present in manifest")
        else:
            expected = entry["total_spots"]
            matches = None
            if expected.isdigit() and r1["read_count"] is not None:
                matches = r1["read_count"] == int(expected)
                if not matches and enforce_manifest:
                    findings.append(
                        f"R1 read count {r1['read_count']} != manifest "
                        f"total_spots {expected}"
                    )
            manifest_info = {
                "present_in_manifest": True,
                "dataset_id": entry["dataset_id"],
                "expected_total_spots": int(expected) if expected.isdigit() else None,
                "matches_r1_read_count": matches,
                "enforced": enforce_manifest,
            }

    malformed = r1["malformed_records"] + r2["malformed_records"]
    if malformed:
        findings.append(f"{malformed} record(s) with non-'@' header line")

    status = "PASS" if not findings else "FAIL"
    return {
        "accession": accession,
        "status": status,
        "findings": findings,
        "r1": r1,
        "r2": r2,
        "pair": {
            "read_count_equal": read_count_equal,
            "id_sync": id_sync,
            "singleton": singleton,
        },
        "manifest": manifest_info,
    }


def build_report(fastq_dir: str, accessions: list[str], manifest_path: Optional[str],
                 max_examples: int, enforce_manifest: bool) -> dict:
    manifest = load_manifest(manifest_path) if manifest_path else None
    entries = [
        validate_accession(acc, fastq_dir, manifest, max_examples, enforce_manifest)
        for acc in accessions
    ]
    passed = sum(1 for e in entries if e["status"] == "PASS")
    return {
        "schema_version": SCHEMA_VERSION,
        "tool": TOOL_NAME,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "fastq_dir": os.path.abspath(fastq_dir),
        "manifest": os.path.abspath(manifest_path) if manifest_path else None,
        "accessions": entries,
        "summary": {
            "total": len(entries),
            "passed": passed,
            "failed": len(entries) - passed,
            "overall_status": "PASS" if passed == len(entries) and entries else "FAIL",
        },
    }


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=TOOL_NAME,
        description="Validate F26 paired FASTQ gzip archives and emit a JSON report.",
    )
    parser.add_argument(
        "--fastq-dir", required=True,
        help="Directory containing <accession>_1.fastq.gz / <accession>_2.fastq.gz",
    )
    parser.add_argument(
        "--accessions", nargs="+", default=None,
        help="Accessions to validate (default: auto-discover *_1.fastq.gz pairs)",
    )
    parser.add_argument(
        "--manifest", default=None,
        help="F26_DATASET_MANIFEST.tsv for optional spot-count cross-check",
    )
    parser.add_argument(
        "--enforce-manifest-spots", action="store_true",
        help="Treat R1 read count != manifest total_spots as a failure",
    )
    parser.add_argument(
        "--max-examples", type=int, default=10,
        help="Maximum read-ID mismatch examples recorded per accession (default: 10)",
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Write JSON report to this path (default: stdout)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    if not os.path.isdir(args.fastq_dir):
        print(f"ERROR: fastq directory not found: {args.fastq_dir}", file=sys.stderr)
        return 2

    accessions = args.accessions or find_accessions(args.fastq_dir)
    if not accessions:
        print(f"ERROR: no *_1.fastq.gz pairs found in {args.fastq_dir}", file=sys.stderr)
        return 2

    try:
        report = build_report(
            args.fastq_dir, accessions, args.manifest,
            args.max_examples, args.enforce_manifest_spots,
        )
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    payload = json.dumps(report, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(payload + "\n")
        print(f"report written: {args.output}", file=sys.stderr)
    else:
        print(payload)

    return 0 if report["summary"]["overall_status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
