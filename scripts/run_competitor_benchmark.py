#!/usr/bin/env python3
"""Manifest-driven F26 competitor benchmark harness.

Protocol: F26_COMPETITOR_BENCHMARK_CONTRACT.md
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

SCHEMA_VERSION = "f26-benchmark-1"

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "F26_DATASET_MANIFEST.tsv"
DEFAULT_DATA_ROOT = Path(
    "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data"
)

SPADES_BIN = Path("/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/spades.py")
UNICYCLER_BIN = Path("/hpcfs/fhome/caizhh/.conda/envs/BactGenome/bin/unicycler")
SKESA_BIN = Path("/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/skesa")
QUAST_BIN = Path("/hpcfs/fhome/caizhh/.local/bin/quast.py")
BACTGENOME_BIN = SPADES_BIN.parent
LINEABACT_BIN = SKESA_BIN.parent
GNU_TIME = Path("/usr/bin/time")

TOOL_ORDER = ("fastbact", "spades", "skesa", "unicycler")

# Partition qfree advertises MaxTime=03:00:00. Defaults cover one dataset,
# one tool, and one replicate, and leave 1200 s for read hashing, version
# probes, and teardown (8100 + 1500 = 9600). A real run whose worst-case
# product exceeds QFREE_MAX_WALL_S is rejected before any assembler starts.
QFREE_MAX_WALL_S = 3 * 60 * 60
DEFAULT_TOOL_TIMEOUT_S = 2 * 60 * 60 + 15 * 60  # 8100
DEFAULT_QUAST_TIMEOUT_S = 25 * 60  # 1500
MIN_DETERMINISM_REPLICATES = 3

# Replicates pool only when the experiment identity matches. Missing fields
# compare as empty so a fixture that omits them still groups with itself.
DETERMINISM_IDENTITY_FIELDS = (
    "dataset_id",
    "tool",
    "r1",
    "r2",
    "reference",
    "reference_sha256",
    "r1_sha256",
    "r2_sha256",
    "binary_sha256",
    "git_head",
    "threads",
    "memory_gb",
)
VERSION_NEEDLES = {
    "spades": "v4.2.0",
    "skesa": "2.5.1",
    "unicycler": "v0.5.1",
    "quast": "v5.2.0",
}

REQUIRED_MANIFEST_COLUMNS = (
    "dataset_id",
    "sra_accession",
    "bioproject",
    "reference_accession",
    "reference_fasta_path",
    "reference_raw_sha256",
)

# NGA50 may be null when QUAST cannot cover half the reference.
REQUIRED_QUAST_FIELDS = (
    "contigs_ge_200bp",
    "n50",
    "l50",
    "largest_contig",
    "total_length",
    "gc_pct",
    "misassemblies",
    "local_misassemblies",
    "mismatches_per_100kbp",
    "indels_per_100kbp",
    "genome_fraction_pct",
    "duplication_ratio",
)

METRIC_FIELDS = REQUIRED_QUAST_FIELDS[:2] + ("nga50",) + REQUIRED_QUAST_FIELDS[2:]

QUAST_KEY_MAP = {
    "contigs_ge_200bp_explicit": "# contigs (>= 200 bp)",
    "contigs_min_contig": "# contigs",
    "n50": "N50",
    "nga50": "NGA50",
    "l50": "L50",
    "largest_contig": "Largest contig",
    "total_length": "Total length",
    "gc_pct": "GC (%)",
    "misassemblies": "# misassemblies",
    "local_misassemblies": "# local misassemblies",
    "mismatches_per_100kbp": "# mismatches per 100 kbp",
    "indels_per_100kbp": "# indels per 100 kbp",
    "genome_fraction_pct": "Genome fraction (%)",
    "duplication_ratio": "Duplication ratio",
}

SUMMARY_COLUMNS = (
    "dataset_id",
    "tool",
    "replicate",
    "status",
    "exit_code",
    "quast_exit_code",
    "wall_clock_s",
    "cpu_time_s",
    "peak_rss_mb",
    "peak_rss_kib",
    "assembly_sha256",
    *METRIC_FIELDS,
    "stderr_log",
    "command",
)


def tool_env(threads: int) -> dict[str, str]:
    env = os.environ.copy()
    prefix = f"{BACTGENOME_BIN}:{LINEABACT_BIN}"
    env["PATH"] = prefix + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = str(threads)
    env["RAYON_NUM_THREADS"] = str(threads)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["LC_ALL"] = "C"
    env["TZ"] = "UTC"
    env.setdefault("RUST_LOG", "info")
    return env


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")
        missing = [name for name in REQUIRED_MANIFEST_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(f"Manifest missing columns {missing}: {path}")
        rows = []
        for row in reader:
            if not row.get("dataset_id"):
                continue
            rows.append({key: (value or "").strip() for key, value in row.items() if key})
        return rows


def select_dataset(rows: Sequence[Mapping[str, str]], dataset_id: str) -> dict[str, str]:
    matches = [row for row in rows if row["dataset_id"] == dataset_id]
    if not matches:
        known = ", ".join(row["dataset_id"] for row in rows)
        raise ValueError(f"Unknown dataset {dataset_id!r}. Known: {known}")
    if len(matches) > 1:
        raise ValueError(f"Dataset {dataset_id!r} appears {len(matches)} times")
    return dict(matches[0])


def _existing_file(candidates: Sequence[Path]) -> Optional[Path]:
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def resolve_reference(row: Mapping[str, str], data_root: Path, repo_root: Path) -> Path:
    raw = Path(row["reference_fasta_path"])
    if raw.is_absolute():
        return raw
    stripped = raw
    if raw.parts and raw.parts[0] == "data":
        stripped = Path(*raw.parts[1:])
    candidates = [
        data_root / stripped,
        data_root / raw,
        repo_root / raw,
        data_root.parent / raw,
    ]
    found = _existing_file(candidates)
    return found if found is not None else candidates[0]


def resolve_reads(row: Mapping[str, str], data_root: Path) -> tuple[Path, Path]:
    explicit_r1 = row.get("fastq_r1", "")
    explicit_r2 = row.get("fastq_r2", "")
    if explicit_r1 and explicit_r2:
        return Path(explicit_r1), Path(explicit_r2)

    acc = row["sra_accession"]
    historical = (
        data_root / "prjna1088182" / "fastq" / acc / f"{acc}_1.fastq.gz",
        data_root / "prjna1088182" / "fastq" / acc / f"{acc}_2.fastq.gz",
    )
    cross = (
        data_root / "f26_cross_species" / "fastq" / f"{acc}_1.fastq.gz",
        data_root / "f26_cross_species" / "fastq" / f"{acc}_2.fastq.gz",
    )
    pairs = [historical, cross]
    if row.get("bioproject") != "PRJNA1088182":
        pairs = [cross, historical]
    for r1, r2 in pairs:
        if r1.is_file() and r2.is_file():
            return r1, r2
    return pairs[0]


def resolve_fastbact_bin(explicit: Optional[str]) -> Path:
    if explicit:
        return Path(explicit)
    candidates = [
        REPO_ROOT / "target" / "release" / "fastBact",
        REPO_ROOT / "target" / "release" / "fastbact",
    ]
    found = _existing_file(candidates)
    return found if found is not None else candidates[0]


def git_head(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return ""
    if result.returncode != 0:
        return ""
    return result.stdout.strip()


def parse_elapsed(text: str) -> Optional[float]:
    token = text.strip()
    if not token:
        return None
    parts = token.split(":")
    try:
        if len(parts) == 1:
            return float(parts[0])
        if len(parts) == 2:
            minutes, seconds = parts
            return int(minutes) * 60 + float(seconds)
        if len(parts) == 3:
            hours, minutes, seconds = parts
            return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    except ValueError:
        return None
    return None


def parse_gnu_time(text: str) -> dict[str, Optional[float]]:
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if ": " not in line:
            continue
        key, value = line.split(": ", 1)
        fields[key.strip()] = value.strip()
    user = _float_or_none(fields.get("User time (seconds)"))
    system = _float_or_none(fields.get("System time (seconds)"))
    cpu = None if user is None or system is None else user + system
    rss_kib = _float_or_none(fields.get("Maximum resident set size (kbytes)"))
    elapsed_key = next((key for key in fields if key.startswith("Elapsed (wall clock) time")), None)
    wall = parse_elapsed(fields[elapsed_key]) if elapsed_key else None
    return {
        "wall_clock_s": wall,
        "cpu_time_s": cpu,
        "peak_rss_kib": rss_kib,
        "peak_rss_mb": None if rss_kib is None else rss_kib / 1024.0,
        "resource_parse_ok": wall is not None and cpu is not None and rss_kib is not None,
    }


def _float_or_none(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_quast_table(text: str) -> dict[str, str]:
    table: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        if key in {"Assembly", ""}:
            continue
        table[key] = parts[1].strip()
    return table


def _metric_number(raw: Optional[str]) -> Any:
    if raw is None:
        return None
    if raw in {"", "-"}:
        return None
    try:
        if any(mark in raw for mark in (".", "e", "E")):
            return float(raw)
        return int(raw)
    except ValueError:
        return raw


def extract_metrics(table: Mapping[str, str]) -> tuple[dict[str, Any], list[str]]:
    notes: list[str] = []
    metrics: dict[str, Any] = {}
    explicit = _metric_number(table.get(QUAST_KEY_MAP["contigs_ge_200bp_explicit"]))
    filtered = _metric_number(table.get(QUAST_KEY_MAP["contigs_min_contig"]))
    if explicit is not None:
        metrics["contigs_ge_200bp"] = explicit
        if filtered is not None and explicit != filtered:
            notes.append(
                f"# contigs (>= 200 bp)={explicit} disagrees with # contigs={filtered}"
            )
    else:
        metrics["contigs_ge_200bp"] = filtered
        if filtered is None:
            notes.append("missing contig count")
    for field, key in QUAST_KEY_MAP.items():
        if field in {"contigs_ge_200bp_explicit", "contigs_min_contig"}:
            continue
        if key not in table:
            metrics[field] = None
            if field != "nga50":
                notes.append(f"missing {key}")
            continue
        metrics[field] = _metric_number(table[key])
        if field != "nga50" and metrics[field] is None and table[key] not in {"-"}:
            notes.append(f"unparsed {key}={table[key]!r}")
    return metrics, notes


def build_tool_command(
    tool: str,
    r1: Path,
    r2: Path,
    tool_out: Path,
    threads: int,
    memory_gb: int,
    fastbact_bin: Path,
) -> list[str]:
    if tool == "fastbact":
        return [
            str(fastbact_bin),
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(tool_out),
            "-t",
            str(threads),
            "--min-len",
            "500",
            "--kmers",
            "auto",
            "--depth",
            "120",
            "--dump-stages",
            str(tool_out / "stages"),
        ]
    if tool == "spades":
        return [
            str(SPADES_BIN),
            "--isolate",
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(tool_out),
            "-t",
            str(threads),
            "-m",
            str(memory_gb),
        ]
    if tool == "skesa":
        return [
            str(SKESA_BIN),
            "--reads",
            f"{r1},{r2}",
            "--cores",
            str(threads),
            "--memory",
            str(memory_gb),
            "--min_contig",
            "200",
            "--contigs_out",
            str(tool_out / "contigs.fa"),
        ]
    if tool == "unicycler":
        return [
            str(UNICYCLER_BIN),
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(tool_out),
            "-t",
            str(threads),
            "--mode",
            "normal",
            "--min_fasta_length",
            "200",
            "--spades_path",
            str(SPADES_BIN),
            "--spades_options",
            f"-m {memory_gb}",
        ]
    raise ValueError(f"Unknown tool {tool!r}")


def build_quast_command(
    reference: Path,
    assembly: Path,
    quast_dir: Path,
    threads: int,
    label: str,
) -> list[str]:
    return [
        str(QUAST_BIN),
        "-R",
        str(reference),
        "-o",
        str(quast_dir),
        "-t",
        str(threads),
        "--min-contig",
        "200",
        "--contig-thresholds",
        "0,200,1000,5000,10000,25000,50000",
        "--labels",
        label,
        str(assembly),
    ]


def tool_binary(tool: str, fastbact_bin: Path) -> Path:
    return {
        "fastbact": fastbact_bin,
        "spades": SPADES_BIN,
        "skesa": SKESA_BIN,
        "unicycler": UNICYCLER_BIN,
    }[tool]


def command_version_text(binary: Path, tool: str) -> tuple[int, str]:
    flag = "--version"
    try:
        result = subprocess.run(
            [str(binary), flag],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, str(exc)
    text = (result.stdout or "") + (result.stderr or "")
    return result.returncode, text.strip()


def version_accepted(tool: str, text: str) -> bool:
    needle = VERSION_NEEDLES.get(tool)
    if needle is None:
        return True
    return needle in text


def classify_exit(exit_code: int, timed_out: bool) -> str:
    if timed_out or exit_code == 124:
        return "TIMEOUT"
    if exit_code == 137:
        return "OOM_KILLED"
    if exit_code == 0:
        return "PASS"
    return "FAIL"


def run_timed(
    argv: Sequence[str],
    stdout_path: Path,
    stderr_path: Path,
    time_path: Path,
    timeout_s: int,
    env: Mapping[str, str],
) -> dict[str, Any]:
    """Run argv under GNU time. Timeout kills the child session only."""
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    timed_out = False
    start = time.monotonic()
    command = [str(GNU_TIME), "-v", "-o", str(time_path), "--", *argv]
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        process = subprocess.Popen(
            command,
            stdout=stdout,
            stderr=stderr,
            env=dict(env),
            start_new_session=True,
        )
        try:
            exit_code = process.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill_session(process)
            exit_code = 124
            stderr.write(
                f"\n[harness] TIMEOUT after {timeout_s}s: process group terminated\n"
            )
    wall_fallback = time.monotonic() - start
    time_text = time_path.read_text(encoding="utf-8", errors="replace") if time_path.exists() else ""
    resources = parse_gnu_time(time_text)
    if resources["wall_clock_s"] is None:
        resources["wall_clock_s"] = wall_fallback
        resources["resource_parse_ok"] = False
    resources.update(
        {
            "exit_code": exit_code,
            "timed_out": timed_out,
            "status": classify_exit(exit_code, timed_out),
            "argv": list(argv),
        }
    )
    return resources


def _kill_session(process: subprocess.Popen[Any]) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=30)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    process.wait(timeout=30)


def locate_assembly(tool_out: Path) -> Optional[Path]:
    for name in ("contigs.fa", "contigs.fasta", "assembly.fasta"):
        candidate = tool_out / name
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return None


def blank_metrics() -> dict[str, Any]:
    return {field: None for field in METRIC_FIELDS}


def run_record(
    *,
    dataset_id: str,
    tool: str,
    replicate: int,
    status: str,
    command: Sequence[str],
    stderr_log: Path,
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": dataset_id,
        "tool": tool,
        "replicate": replicate,
        "status": status,
        "exit_code": None,
        "quast_exit_code": None,
        "wall_clock_s": None,
        "cpu_time_s": None,
        "peak_rss_mb": None,
        "peak_rss_kib": None,
        "quast_wall_clock_s": None,
        "quast_cpu_time_s": None,
        "quast_peak_rss_mb": None,
        "assembly_sha256": None,
        "command": list(command),
        "stderr_log": str(stderr_log),
        "metrics": blank_metrics(),
        "notes": [],
    }
    if extra:
        record.update(extra)
    return record


def execute_tool(
    *,
    dataset_id: str,
    tool: str,
    replicate: int,
    r1: Path,
    r2: Path,
    reference: Path,
    run_dir: Path,
    threads: int,
    memory_gb: int,
    timeout_s: int,
    quast_timeout_s: int,
    fastbact_bin: Path,
    dry_run: bool,
    force: bool,
) -> dict[str, Any]:
    tool_out = run_dir / "tool_out"
    stderr_log = run_dir / "stderr.log"
    command = build_tool_command(tool, r1, r2, tool_out, threads, memory_gb, fastbact_bin)
    label = f"{dataset_id}-{tool}-rep{replicate}"

    if run_dir.exists() and not dry_run and not force:
        record = run_record(
            dataset_id=dataset_id,
            tool=tool,
            replicate=replicate,
            status="EXISTS",
            command=command,
            stderr_log=stderr_log,
        )
        record["notes"] = [f"refusing to overwrite {run_dir}; pass --force"]
        return record

    if dry_run:
        quast_command = build_quast_command(
            reference, run_dir / "contigs.fa", run_dir / "quast", threads, label
        )
        return run_record(
            dataset_id=dataset_id,
            tool=tool,
            replicate=replicate,
            status="DRY_RUN",
            command=command,
            stderr_log=stderr_log,
            extra={"quast_command": quast_command, "r1": str(r1), "r2": str(r2), "reference": str(reference)},
        )

    if run_dir.exists() and force:
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")

    binary = tool_binary(tool, fastbact_bin)
    if not binary.is_file():
        stderr_log.write_text(f"missing binary: {binary}\n", encoding="utf-8")
        record = run_record(
            dataset_id=dataset_id,
            tool=tool,
            replicate=replicate,
            status="MISSING_BINARY",
            command=command,
            stderr_log=stderr_log,
        )
        _write_json(run_dir / "metrics.json", record)
        return record

    if tool != "fastbact":
        version_code, version_text = command_version_text(binary, tool)
        (run_dir / "version.txt").write_text(version_text + "\n", encoding="utf-8")
        if version_code != 0 or not version_accepted(tool, version_text):
            record = run_record(
                dataset_id=dataset_id,
                tool=tool,
                replicate=replicate,
                status="VERSION_MISMATCH",
                command=command,
                stderr_log=stderr_log,
                extra={"tool_version_text": version_text, "version_exit_code": version_code},
            )
            _write_json(run_dir / "metrics.json", record)
            return record
    else:
        version_code, version_text = command_version_text(binary, tool)
        (run_dir / "version.txt").write_text(version_text + "\n", encoding="utf-8")
        if binary.is_file():
            (run_dir / "binary_sha256.txt").write_text(sha256_file(binary) + "\n", encoding="utf-8")

    try:
        timed = run_timed(
            command,
            run_dir / "stdout.log",
            stderr_log,
            run_dir / "time.txt",
            timeout_s,
            tool_env(threads),
        )
    except Exception as exc:  # noqa: BLE001 — one tool must not abort the sweep
        stderr_log.write_text(traceback.format_exc(), encoding="utf-8")
        record = run_record(
            dataset_id=dataset_id,
            tool=tool,
            replicate=replicate,
            status="FAIL",
            command=command,
            stderr_log=stderr_log,
            extra={"notes": [f"harness exception: {exc}"]},
        )
        _write_json(run_dir / "metrics.json", record)
        return record

    record = run_record(
        dataset_id=dataset_id,
        tool=tool,
        replicate=replicate,
        status=str(timed["status"]),
        command=command,
        stderr_log=stderr_log,
        extra={
            "exit_code": timed["exit_code"],
            "wall_clock_s": timed["wall_clock_s"],
            "cpu_time_s": timed["cpu_time_s"],
            "peak_rss_mb": timed["peak_rss_mb"],
            "peak_rss_kib": timed["peak_rss_kib"],
            "resource_parse_ok": timed["resource_parse_ok"],
            "tool_version_text": version_text,
        },
    )
    if record["status"] != "PASS":
        _write_json(run_dir / "metrics.json", record)
        return record

    produced = locate_assembly(tool_out)
    if produced is None:
        record["status"] = "MISSING_ASSEMBLY"
        record["notes"].append("assembler exit 0 but no contig FASTA")
        _write_json(run_dir / "metrics.json", record)
        return record

    canonical = run_dir / "contigs.fa"
    shutil.copy2(produced, canonical)
    record["assembly_sha256"] = sha256_file(canonical)
    record["assembly_source"] = str(produced)

    quast_dir = run_dir / "quast"
    quast_command = build_quast_command(reference, canonical, quast_dir, threads, label)
    record["quast_command"] = quast_command
    try:
        quast_timed = run_timed(
            quast_command,
            run_dir / "quast.stdout.log",
            run_dir / "quast.stderr.log",
            run_dir / "quast.time.txt",
            quast_timeout_s,
            tool_env(threads),
        )
    except Exception as exc:  # noqa: BLE001
        (run_dir / "quast.stderr.log").write_text(traceback.format_exc(), encoding="utf-8")
        record["status"] = "QUAST_FAIL"
        record["notes"].append(f"quast harness exception: {exc}")
        _write_json(run_dir / "metrics.json", record)
        return record

    record["quast_exit_code"] = quast_timed["exit_code"]
    record["quast_wall_clock_s"] = quast_timed["wall_clock_s"]
    record["quast_cpu_time_s"] = quast_timed["cpu_time_s"]
    record["quast_peak_rss_mb"] = quast_timed["peak_rss_mb"]
    report = quast_dir / "report.tsv"
    if quast_timed["status"] != "PASS" or not report.is_file():
        record["status"] = "TIMEOUT" if quast_timed["status"] == "TIMEOUT" else "QUAST_FAIL"
        record["notes"].append(f"quast status {quast_timed['status']} exit {quast_timed['exit_code']}")
        _write_json(run_dir / "metrics.json", record)
        return record

    table = parse_quast_table(report.read_text(encoding="utf-8", errors="replace"))
    metrics, notes = extract_metrics(table)
    record["metrics"] = metrics
    record["notes"].extend(notes)
    if notes:
        record["status"] = "METRICS_INCOMPLETE"
    _write_json(run_dir / "metrics.json", record)
    return record


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def flatten_row(record: Mapping[str, Any]) -> dict[str, Any]:
    metrics = record.get("metrics") or {}
    flat: dict[str, Any] = {}
    for column in SUMMARY_COLUMNS:
        if column in METRIC_FIELDS:
            flat[column] = metrics.get(column)
        elif column == "command":
            command = record.get("command") or []
            flat[column] = " ".join(str(part) for part in command)
        else:
            flat[column] = record.get(column)
    return flat


def _identity_token(record: Mapping[str, Any], field: str) -> str:
    value = record.get(field)
    if value is None:
        return ""
    return str(value)


def determinism_group_key(record: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(_identity_token(record, field) for field in DETERMINISM_IDENTITY_FIELDS)


def determinism_identity(record: Mapping[str, Any]) -> dict[str, str]:
    return {
        field: _identity_token(record, field)
        for field in DETERMINISM_IDENTITY_FIELDS
        if field not in {"dataset_id", "tool"}
    }


def scheduled_wall_seconds(
    n_tools: int,
    n_replicates: int,
    timeout_s: int,
    quast_timeout_s: int,
) -> int:
    return n_tools * n_replicates * (timeout_s + quast_timeout_s)


def qfree_schedule_error(
    n_tools: int,
    n_replicates: int,
    timeout_s: int,
    quast_timeout_s: int,
) -> Optional[str]:
    """Return an error when the worst-case schedule cannot finish inside qfree."""
    budget = scheduled_wall_seconds(n_tools, n_replicates, timeout_s, quast_timeout_s)
    if budget <= QFREE_MAX_WALL_S:
        return None
    return (
        f"scheduled worst-case wall {budget}s exceeds qfree MaxTime "
        f"{QFREE_MAX_WALL_S}s ({n_tools} tool(s) x {n_replicates} replicate(s) x "
        f"({timeout_s}s assembler + {quast_timeout_s}s QUAST)). "
        "Submit one dataset and one tool per Slurm job. "
        "Certify N=3 with three jobs that share --outdir and pass "
        "--replicate 1, --replicate 2, and --replicate 3."
    )


def replicate_indices(args: argparse.Namespace) -> list[int]:
    if getattr(args, "replicate", None) is not None:
        return [int(args.replicate)]
    return list(range(1, int(args.replicates) + 1))


def stamp_run_identity(
    record: dict[str, Any],
    *,
    r1: Path,
    r2: Path,
    reference: Path,
    reference_sha256: str,
    r1_sha256: str,
    r2_sha256: str,
    binary_sha256: str,
    git_head_value: str,
    threads: int,
    memory_gb: int,
) -> None:
    record["r1"] = str(r1)
    record["r2"] = str(r2)
    record["reference"] = str(reference)
    record["reference_sha256"] = reference_sha256
    record["r1_sha256"] = r1_sha256
    record["r2_sha256"] = r2_sha256
    record["binary_sha256"] = binary_sha256
    record["git_head"] = git_head_value
    record["threads"] = threads
    record["memory_gb"] = memory_gb


def determinism_report(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault(determinism_group_key(record), []).append(record)
    reports = []
    for key, group in sorted(grouped.items()):
        dataset_id, tool = key[0], key[1]
        identity = determinism_identity(group[0])
        if len(group) < MIN_DETERMINISM_REPLICATES:
            reports.append(
                {
                    "dataset_id": dataset_id,
                    "tool": tool,
                    "replicates": len(group),
                    "checked": False,
                    "pass": None,
                    "reason": (
                        f"not enough replicates: {len(group)} < {MIN_DETERMINISM_REPLICATES}"
                    ),
                    "identity": identity,
                }
            )
            continue
        statuses = [str(item["status"]) for item in group]
        hashes = [item.get("assembly_sha256") for item in group]
        metric_rows = [item.get("metrics") or {} for item in group]
        mismatches = []
        if any(status != "PASS" for status in statuses):
            mismatches.append("non-PASS replicate")
        if len(set(hashes)) != 1 or hashes[0] in {None, ""}:
            mismatches.append("assembly_sha256 differs or is null")
        variance_max = 0.0
        for field in METRIC_FIELDS:
            values = [row.get(field) for row in metric_rows]
            if any(value != values[0] for value in values):
                mismatches.append(field)
            numeric = [value for value in values if isinstance(value, (int, float))]
            if len(numeric) == len(values) and numeric:
                center = sum(numeric) / len(numeric)
                variance = sum((value - center) ** 2 for value in numeric) / len(numeric)
                variance_max = max(variance_max, variance)
        reports.append(
            {
                "dataset_id": dataset_id,
                "tool": tool,
                "replicates": len(group),
                "checked": True,
                "pass": not mismatches,
                "variance_max": variance_max,
                "assembly_sha256": hashes[0] if len(set(hashes)) == 1 else None,
                "mismatches": mismatches,
                "identity": identity,
            }
        )
    return reports


def load_existing_summary(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        backup = path.with_suffix(".json.corrupt")
        shutil.copy2(path, backup)
        print(f"[harness] preserved unreadable summary at {backup}", file=sys.stderr)
        return []
    runs = payload.get("runs") if isinstance(payload, dict) else None
    if not isinstance(runs, list):
        return []
    return [item for item in runs if isinstance(item, dict)]


def merge_runs(existing: Sequence[Mapping[str, Any]], new_runs: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    def key(record: Mapping[str, Any]) -> tuple[str, str, int]:
        return (str(record.get("dataset_id")), str(record.get("tool")), int(record.get("replicate") or 0))

    merged = {key(record): dict(record) for record in existing}
    for record in new_runs:
        merged[key(record)] = dict(record)
    return list(merged.values())


def write_summaries(outdir: Path, runs: Sequence[Mapping[str, Any]], context: Mapping[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    reports = determinism_report(runs)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "context": dict(context),
        "runs": list(runs),
        "determinism": reports,
    }
    _write_json(outdir / "summary.json", payload)
    with (outdir / "summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SUMMARY_COLUMNS), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in runs:
            writer.writerow(flatten_row(record))


def input_failure(
    *,
    dataset_id: str,
    tools: Sequence[str],
    replicate_ids: Sequence[int],
    status: str,
    note: str,
    r1: Path,
    r2: Path,
    reference: Path,
    fastbact_bin: Path,
    threads: int,
    memory_gb: int,
) -> list[dict[str, Any]]:
    records = []
    for tool in tools:
        for replicate in replicate_ids:
            command = build_tool_command(
                tool,
                r1,
                r2,
                Path(f"{{outdir}}/{dataset_id}/{tool}/rep_{replicate}/tool_out"),
                threads,
                memory_gb,
                fastbact_bin,
            )
            records.append(
                run_record(
                    dataset_id=dataset_id,
                    tool=tool,
                    replicate=replicate,
                    status=status,
                    command=command,
                    stderr_log=Path(f"{{outdir}}/{dataset_id}/{tool}/rep_{replicate}/stderr.log"),
                    extra={"notes": [note], "r1": str(r1), "r2": str(r2), "reference": str(reference)},
                )
            )
    return records


def preflight_inputs(
    row: Mapping[str, str],
    r1: Path,
    r2: Path,
    reference: Path,
    verify_hash: bool,
) -> tuple[Optional[str], str]:
    if not r1.is_file() or not r2.is_file():
        return "MISSING_INPUT", f"read pair missing: {r1} , {r2}"
    if not reference.is_file():
        return "MISSING_REFERENCE", f"reference missing: {reference}"
    if verify_hash:
        expected = row["reference_raw_sha256"].lower()
        observed = sha256_file(reference)
        if observed != expected:
            return (
                "REFERENCE_HASH_MISMATCH",
                f"reference sha256 {observed} != manifest {expected}",
            )
    return None, ""


def evaluate_dataset(args: argparse.Namespace, row: Mapping[str, str]) -> int:
    data_root = Path(args.data_root)
    r1, r2 = resolve_reads(row, data_root)
    reference = resolve_reference(row, data_root, REPO_ROOT)
    fastbact_bin = resolve_fastbact_bin(args.fastbact_bin)
    tools = list(TOOL_ORDER) if args.tool == "all" else [args.tool]
    outdir = Path(args.outdir)
    context = {
        "dataset_id": row["dataset_id"],
        "sra_accession": row["sra_accession"],
        "reference_accession": row["reference_accession"],
        "reference": str(reference),
        "r1": str(r1),
        "r2": str(r2),
        "threads": args.threads,
        "memory_gb": args.memory_gb,
        "git_head": git_head(REPO_ROOT),
        "fastbact_bin": str(fastbact_bin),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "hostname": os.uname().nodename,
        "dry_run": bool(args.dry_run),
    }
    if args.threads != 8 or args.memory_gb != 14:
        context["limit_deviation"] = (
            f"protocol is 8 threads and 14 GB; this invocation uses "
            f"{args.threads} threads and {args.memory_gb} GB"
        )
        print(f"[harness] {context['limit_deviation']}", file=sys.stderr)

    replicate_ids = replicate_indices(args)
    schedule_error = qfree_schedule_error(
        len(tools),
        len(replicate_ids),
        args.timeout_seconds,
        args.quast_timeout_seconds,
    )
    if schedule_error:
        print(f"[harness] {schedule_error}", file=sys.stderr)
        if not args.dry_run:
            return 2

    failure, note = preflight_inputs(row, r1, r2, reference, not args.no_verify_reference_hash and not args.dry_run)
    if args.dry_run:
        failure, note = preflight_inputs(row, r1, r2, reference, not args.no_verify_reference_hash)

    r1_sha256 = ""
    r2_sha256 = ""
    if not failure and not args.dry_run:
        print("[harness] hashing reads for determinism identity", flush=True)
        r1_sha256 = sha256_file(r1)
        r2_sha256 = sha256_file(r2)
    binary_hashes: dict[str, str] = {}
    for tool in tools:
        binary = tool_binary(tool, fastbact_bin)
        binary_hashes[tool] = sha256_file(binary) if binary.is_file() else ""

    new_runs: list[dict[str, Any]]
    if failure and not args.dry_run:
        print(f"[harness] {failure}: {note}", file=sys.stderr)
        new_runs = input_failure(
            dataset_id=row["dataset_id"],
            tools=tools,
            replicate_ids=replicate_ids,
            status=failure,
            note=note,
            r1=r1,
            r2=r2,
            reference=reference,
            fastbact_bin=fastbact_bin,
            threads=args.threads,
            memory_gb=args.memory_gb,
        )
    else:
        new_runs = []
        for tool in tools:
            for replicate in replicate_ids:
                run_dir = outdir / row["dataset_id"] / tool / f"rep_{replicate}"
                print(f"[harness] {row['dataset_id']} {tool} replicate {replicate} -> {run_dir}", flush=True)
                try:
                    record = execute_tool(
                        dataset_id=row["dataset_id"],
                        tool=tool,
                        replicate=replicate,
                        r1=r1,
                        r2=r2,
                        reference=reference,
                        run_dir=run_dir,
                        threads=args.threads,
                        memory_gb=args.memory_gb,
                        timeout_s=args.timeout_seconds,
                        quast_timeout_s=args.quast_timeout_seconds,
                        fastbact_bin=fastbact_bin,
                        dry_run=bool(args.dry_run),
                        force=bool(args.force),
                    )
                except Exception as exc:  # noqa: BLE001
                    print(f"[harness] isolated failure on {tool}: {exc}", file=sys.stderr)
                    record = run_record(
                        dataset_id=row["dataset_id"],
                        tool=tool,
                        replicate=replicate,
                        status="FAIL",
                        command=[],
                        stderr_log=run_dir / "stderr.log",
                        extra={"notes": [traceback.format_exc()]},
                    )
                if failure and args.dry_run:
                    record["status"] = failure
                    record.setdefault("notes", []).append(note)
                new_runs.append(record)
                print(f"[harness] status={record['status']} exit={record['exit_code']}", flush=True)

    reference_sha256 = row["reference_raw_sha256"].lower()
    for record in new_runs:
        stamp_run_identity(
            record,
            r1=r1,
            r2=r2,
            reference=reference,
            reference_sha256=reference_sha256,
            r1_sha256=r1_sha256,
            r2_sha256=r2_sha256,
            binary_sha256=binary_hashes[str(record["tool"])],
            git_head_value=str(context["git_head"]),
            threads=args.threads,
            memory_gb=args.memory_gb,
        )

    merged = merge_runs(load_existing_summary(outdir / "summary.json"), new_runs)
    write_summaries(outdir, merged, context)
    print(f"[harness] wrote {outdir / 'summary.json'}", flush=True)
    failed = [record for record in new_runs if record["status"] not in {"PASS", "DRY_RUN"}]
    return 1 if failed else 0


def self_test() -> int:
    failures: list[str] = []

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)
            print(f"FAIL {message}")
        else:
            print(f"ok   {message}")

    elapsed = parse_elapsed("1:02:03.5")
    check(elapsed == 3723.5, f"elapsed h:mm:ss.s -> {elapsed}")
    check(parse_elapsed("3:01.25") == 181.25, "elapsed m:ss.ss")

    time_text = (
        "Command being timed: \"sleep\"\n"
        "User time (seconds): 1.5\n"
        "System time (seconds): 0.25\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:02.00\n"
        "Maximum resident set size (kbytes): 2048\n"
    )
    parsed = parse_gnu_time(time_text)
    check(parsed["cpu_time_s"] == 1.75, "cpu time sums user and system")
    check(parsed["wall_clock_s"] == 2.0, "wall clock parsed")
    check(parsed["peak_rss_mb"] == 2.0, "peak RSS converted to MiB")

    report = (
        "Assembly\tdemo\n"
        "# contigs (>= 0 bp)\t10\n"
        "# contigs (>= 200 bp)\t4\n"
        "# contigs\t4\n"
        "Largest contig\t1000\n"
        "Total length\t2500\n"
        "GC (%)\t50.5\n"
        "N50\t800\n"
        "L50\t2\n"
        "# misassemblies\t1\n"
        "# local misassemblies\t3\n"
        "# mismatches per 100 kbp\t0.85\n"
        "# indels per 100 kbp\t-\n"
        "Genome fraction (%)\t98.115\n"
        "Duplication ratio\t1.001\n"
        "NGA50\t-\n"
    )
    metrics, notes = extract_metrics(parse_quast_table(report))
    check(metrics["contigs_ge_200bp"] == 4, "contigs >= 200")
    check(metrics["n50"] == 800 and metrics["l50"] == 2, "N50 and L50")
    check(metrics["gc_pct"] == 50.5, "GC")
    check(metrics["misassemblies"] == 1 and metrics["local_misassemblies"] == 3, "misassemblies")
    check(metrics["indels_per_100kbp"] is None and metrics["nga50"] is None, "dash is null")
    check(notes == [], f"complete report has no notes, got {notes}")

    disagree, disagree_notes = extract_metrics(
        {"# contigs (>= 200 bp)": "3", "# contigs": "4", "N50": "1", "L50": "1",
         "Largest contig": "1", "Total length": "1", "GC (%)": "1", "# misassemblies": "0",
         "# local misassemblies": "0", "# mismatches per 100 kbp": "0",
         "# indels per 100 kbp": "0", "Genome fraction (%)": "1", "Duplication ratio": "1"}
    )
    check(any("disagrees" in note for note in disagree_notes), "contig row disagreement is flagged")
    check(disagree["contigs_ge_200bp"] == 3, "explicit 200 bp row wins")

    rows = load_manifest(DEFAULT_MANIFEST)
    row = select_dataset(rows, "R1-001")
    r1, r2 = resolve_reads(row, DEFAULT_DATA_ROOT)
    reference = resolve_reference(row, DEFAULT_DATA_ROOT, REPO_ROOT)
    check(r1.is_file() and r2.is_file(), f"R1-001 reads resolve to {r1}")
    check(reference.is_file(), f"R1-001 reference resolves to {reference}")
    check(sha256_file(reference) == row["reference_raw_sha256"], "R1-001 reference sha256 matches manifest")

    cross = select_dataset(rows, "CS-SA01")
    cross_r1, _cross_r2 = resolve_reads(cross, DEFAULT_DATA_ROOT)
    check(
        "f26_cross_species" in str(cross_r1),
        f"CS-SA01 canonical read path uses cross-species layout ({cross_r1})",
    )

    command = build_tool_command(
        "skesa", r1, r2, Path("/tmp/tool_out"), 8, 14, Path("/tmp/fastBact")
    )
    check(command[command.index("--reads") + 1] == f"{r1},{r2}", "SKESA reads are one comma-joined argument")
    check("--min_contig" in command and command[command.index("--min_contig") + 1] == "200", "SKESA min_contig 200")
    spades = build_tool_command("spades", r1, r2, Path("/tmp/tool_out"), 8, 14, Path("/tmp/fastBact"))
    check(spades[1] == "--isolate" and "-m" in spades and spades[spades.index("-m") + 1] == "14", "SPAdes isolate -m 14")
    uni = build_tool_command("unicycler", r1, r2, Path("/tmp/tool_out"), 8, 14, Path("/tmp/fastBact"))
    check("--mode" in uni and uni[uni.index("--spades_options") + 1] == "-m 14", "Unicycler pins SPAdes memory")
    check("-l" not in uni, "Unicycler command has no long-read input")
    fb = build_tool_command("fastbact", r1, r2, Path("/tmp/tool_out"), 8, 14, Path("/tmp/fastBact"))
    check("--enable-auto-subsample" not in fb and "--min-len" in fb, "fastBact keeps F25 defaults without subsampling")
    quast = build_quast_command(reference, Path("/tmp/contigs.fa"), Path("/tmp/quast"), 8, "R1-001-fastbact-rep1")
    check(quast[1] == "-R" and "--min-contig" in quast and "200" in quast, "QUAST -R and min-contig 200")

    tmp = Path("/tmp/f26-harness-self-test")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    ok_run = run_timed(
        ["/bin/sleep", "0.1"],
        tmp / "ok.out",
        tmp / "ok.err",
        tmp / "ok.time",
        30,
        tool_env(1),
    )
    check(ok_run["status"] == "PASS" and ok_run["exit_code"] == 0, f"sleep succeeds ({ok_run['status']})")
    check(ok_run["peak_rss_kib"] is not None and ok_run["wall_clock_s"] is not None, "GNU time fields captured")
    fail_run = run_timed(
        ["/bin/false"],
        tmp / "fail.out",
        tmp / "fail.err",
        tmp / "fail.time",
        30,
        tool_env(1),
    )
    check(fail_run["status"] == "FAIL" and fail_run["exit_code"] != 0, "false is FAIL and does not raise")
    timeout_run = run_timed(
        ["/bin/sleep", "30"],
        tmp / "timeout.out",
        tmp / "timeout.err",
        tmp / "timeout.time",
        1,
        tool_env(1),
    )
    check(timeout_run["status"] == "TIMEOUT" and timeout_run["exit_code"] == 124, "timeout is recorded")
    check("TIMEOUT" in (tmp / "timeout.err").read_text(encoding="utf-8"), "timeout note written to stderr")

    def one_replicate(**overrides: Any) -> dict[str, Any]:
        record: dict[str, Any] = {
            "dataset_id": "R1-001",
            "tool": "fastbact",
            "status": "PASS",
            "assembly_sha256": "abc",
            "r1": "/reads/R1.fastq.gz",
            "r2": "/reads/R2.fastq.gz",
            "reference": "/ref/genome.fasta",
            "reference_sha256": "refhash",
            "r1_sha256": "r1hash",
            "r2_sha256": "r2hash",
            "binary_sha256": "binhash",
            "git_head": "abc123",
            "threads": 8,
            "memory_gb": 14,
            "metrics": {field: 1 for field in METRIC_FIELDS},
        }
        record.update(overrides)
        return record

    det = determinism_report([one_replicate() for _ in range(3)])
    check(
        len(det) == 1
        and det[0]["replicates"] == 3
        and det[0]["checked"] is True
        and det[0]["pass"] is True
        and det[0]["variance_max"] == 0.0,
        "three identical replicates pass determinism",
    )
    short = determinism_report([one_replicate() for _ in range(2)])
    check(
        short[0]["checked"] is False
        and short[0]["pass"] is None
        and short[0]["replicates"] == 2
        and "not enough replicates" in short[0]["reason"],
        "two replicates are not a determinism PASS",
    )
    det_fail = determinism_report(
        [
            one_replicate(tool="spades"),
            one_replicate(tool="spades"),
            one_replicate(tool="spades", assembly_sha256="def"),
        ]
    )
    check(det_fail[0]["checked"] is True and det_fail[0]["pass"] is False, "hash mismatch fails determinism")
    split = determinism_report(
        [
            one_replicate(),
            one_replicate(),
            one_replicate(r1="/reads/other_R1.fastq.gz", r1_sha256="other"),
        ]
    )
    check(
        all(item["pass"] is not True for item in split),
        "different read identity is not pooled into a determinism PASS",
    )
    check(
        scheduled_wall_seconds(1, 1, DEFAULT_TOOL_TIMEOUT_S, DEFAULT_QUAST_TIMEOUT_S) <= QFREE_MAX_WALL_S,
        "default one-tool schedule fits qfree",
    )
    check(
        qfree_schedule_error(len(TOOL_ORDER), 1, DEFAULT_TOOL_TIMEOUT_S, DEFAULT_QUAST_TIMEOUT_S) is not None,
        "default four-tool schedule exceeds qfree",
    )
    check(
        qfree_schedule_error(1, 3, DEFAULT_TOOL_TIMEOUT_S, DEFAULT_QUAST_TIMEOUT_S) is not None,
        "three default-timeout replicates exceed one qfree job",
    )
    shutil.rmtree(tmp)

    if failures:
        print(f"SELF_TEST_FAIL {len(failures)}")
        return 1
    print("SELF_TEST_PASS")
    return 0


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="F26 unified competitor benchmark harness")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--dataset", help="Dataset ID from the manifest")
    parser.add_argument("--tool", choices=[*TOOL_ORDER, "all"], help="Assembler or all")
    parser.add_argument("--outdir", type=Path, help="Directory for logs, QUAST, and summaries")
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument(
        "--replicate",
        type=int,
        default=None,
        help=(
            "Run one 1-based replicate index. A qfree job is one dataset and one tool; "
            "N=3 determinism is three jobs (--replicate 1, 2, and 3) sharing --outdir."
        ),
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=int, default=14)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TOOL_TIMEOUT_S)
    parser.add_argument("--quast-timeout-seconds", type=int, default=DEFAULT_QUAST_TIMEOUT_S)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--fastbact-bin", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-verify-reference-hash", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        return args
    missing = [name for name in ("dataset", "tool", "outdir") if not getattr(args, name)]
    if missing:
        parser.error("required unless --self-test: " + ", ".join("--" + name for name in missing))
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    if args.replicate is not None:
        if args.replicate < 1:
            parser.error("--replicate must be >= 1")
        if args.replicates != 1:
            parser.error("--replicate selects one index; do not combine it with --replicates greater than 1")
    if args.threads < 1 or args.memory_gb < 1:
        parser.error("--threads and --memory-gb must be >= 1")
    if args.timeout_seconds < 1 or args.quast_timeout_seconds < 1:
        parser.error("timeouts must be >= 1")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        return self_test()
    try:
        rows = load_manifest(args.manifest)
        row = select_dataset(rows, args.dataset)
    except (OSError, ValueError) as exc:
        print(f"[harness] {exc}", file=sys.stderr)
        return 2
    return evaluate_dataset(args, row)


if __name__ == "__main__":
    sys.exit(main())
