#!/usr/bin/env python3
"""Collate isolated F26 benchmark runs into one summary.

Each Stage D Slurm job writes its own directory
(benchmarks/f26_results/raw/${dataset}_${tool}_r1/) so concurrent jobs do not
read-modify-write the same summary.json. This script finds every metrics.json
under --raw-dir and rebuilds summary.json and summary.tsv with the harness
schema (merge_runs, flatten_row, write_summaries).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_competitor_benchmark as harness  # noqa: E402


class CollationError(Exception):
    """A metrics file or raw directory cannot be collated."""


def discover_metrics(raw_dir: Path) -> list[Path]:
    if not raw_dir.is_dir():
        raise CollationError(f"raw dir is not a directory: {raw_dir}")
    return sorted(path for path in raw_dir.rglob("metrics.json") if path.is_file())


def load_record(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CollationError(f"unreadable metrics file {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise CollationError(f"metrics file is not an object: {path}")
    missing = [key for key in ("dataset_id", "tool", "replicate") if key not in payload]
    if missing:
        raise CollationError(f"{path} missing {', '.join(missing)}")
    try:
        replicate = int(payload["replicate"])
    except (TypeError, ValueError) as exc:
        raise CollationError(f"{path} replicate is not an integer") from exc
    record = dict(payload)
    record["replicate"] = replicate
    return record


def _tool_rank(tool: str) -> tuple[int, str]:
    try:
        return (harness.TOOL_ORDER.index(tool), tool)
    except ValueError:
        return (len(harness.TOOL_ORDER), tool)


def ordered_records(raw_dir: Path) -> tuple[list[tuple[Path, dict[str, Any]]], list[dict[str, Any]]]:
    loaded = [(path, load_record(path)) for path in discover_metrics(raw_dir)]
    loaded.sort(
        key=lambda item: (
            str(item[1].get("dataset_id")),
            _tool_rank(str(item[1].get("tool"))),
            int(item[1]["replicate"]),
            str(item[0]),
        )
    )
    merged = harness.merge_runs([], [record for _path, record in loaded])
    return loaded, merged


def collation_context(
    raw_dir: Path,
    loaded: Sequence[tuple[Path, Mapping[str, Any]]],
    runs: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    return {
        "producer": "scripts/collate_benchmark_results.py",
        "raw_dir": str(raw_dir),
        "n_metrics_files": len(loaded),
        "n_runs": len(runs),
        "source_metrics": [str(path.relative_to(raw_dir)) for path, _record in loaded],
        "git_head": harness.git_head(harness.REPO_ROOT),
    }


def collate(raw_dir: Path, outdir: Path, dry_run: bool) -> int:
    raw_dir = raw_dir.resolve()
    loaded, runs = ordered_records(raw_dir)
    context = collation_context(raw_dir, loaded, runs)
    if not loaded:
        print(f"[collate] no metrics.json under {raw_dir}", file=sys.stderr, flush=True)
    if dry_run:
        print(f"[collate] dry-run: {len(loaded)} metrics.json -> {len(runs)} runs", flush=True)
        print(f"[collate] would write {outdir / 'summary.json'}", flush=True)
        print(f"[collate] would write {outdir / 'summary.tsv'}", flush=True)
        for path, record in loaded:
            print(
                "[collate] "
                f"{record.get('dataset_id')} {record.get('tool')} "
                f"replicate {record.get('replicate')} <- {path.relative_to(raw_dir)}",
                flush=True,
            )
        return 0
    harness.write_summaries(outdir, runs, context)
    print(f"[collate] wrote {len(runs)} runs to {outdir / 'summary.json'}", flush=True)
    return 0


def _fixture(dataset_id: str, tool: str, replicate: int, **overrides: Any) -> dict[str, Any]:
    metrics = {field: 1 for field in harness.METRIC_FIELDS}
    metrics["nga50"] = None
    record: dict[str, Any] = {
        "schema_version": harness.SCHEMA_VERSION,
        "dataset_id": dataset_id,
        "tool": tool,
        "replicate": replicate,
        "status": "PASS",
        "exit_code": 0,
        "quast_exit_code": 0,
        "wall_clock_s": 1.5,
        "cpu_time_s": 2.0,
        "peak_rss_mb": 10.0,
        "peak_rss_kib": 10240,
        "assembly_sha256": "abc",
        "command": [tool, "--example"],
        "stderr_log": "stderr.log",
        "metrics": metrics,
        "notes": [],
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
    }
    record.update(overrides)
    return record


def _write_metrics(path: Path, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def self_test() -> int:
    failures: list[str] = []

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)
            print(f"FAIL {message}")
        else:
            print(f"ok   {message}")

    with tempfile.TemporaryDirectory(prefix="f26-collate-self-test-") as tmp_name:
        tmp = Path(tmp_name)
        raw = tmp / "raw"
        _write_metrics(
            raw / "CS-BS01_fastbact_r1" / "CS-BS01" / "fastbact" / "rep_1" / "metrics.json",
            _fixture("CS-BS01", "fastbact", 1, assembly_sha256="bs"),
        )
        _write_metrics(
            raw / "CS-SA01_spades_r1" / "CS-SA01" / "spades" / "rep_1" / "metrics.json",
            _fixture("CS-SA01", "spades", 1, assembly_sha256="spades"),
        )
        _write_metrics(
            raw / "aa_old" / "metrics.json",
            _fixture("CS-SA01", "fastbact", 1, assembly_sha256="loser", status="FAIL"),
        )
        _write_metrics(
            raw / "zz_new" / "metrics.json",
            _fixture("CS-SA01", "fastbact", 1, assembly_sha256="winner"),
        )
        (raw / "CS-SA01_spades_r1" / "notes.json").write_text("{}\n", encoding="utf-8")
        (raw / "CS-SA01_spades_r1" / "CS-SA01" / "spades" / "rep_1" / "quast").mkdir(parents=True)
        (raw / "CS-SA01_spades_r1" / "CS-SA01" / "spades" / "rep_1" / "quast" / "report.tsv").write_text(
            "Assembly\tdemo\n",
            encoding="utf-8",
        )

        out = tmp / "out"
        code = main(["--raw-dir", str(raw), "--outdir", str(out)])
        check(code == 0, f"collate exit {code}")
        summary_path = out / "summary.json"
        tsv_path = out / "summary.tsv"
        check(summary_path.is_file() and tsv_path.is_file(), "summary.json and summary.tsv exist")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        check(summary["schema_version"] == harness.SCHEMA_VERSION, "schema_version matches harness")
        check(summary["context"]["n_metrics_files"] == 4, "four metrics.json files discovered")
        check(summary["context"]["n_runs"] == 3, "duplicate dataset/tool/replicate collapsed")
        check(isinstance(summary.get("determinism"), list), "determinism report present")
        runs = summary["runs"]
        check(
            [(item["dataset_id"], item["tool"]) for item in runs]
            == [("CS-BS01", "fastbact"), ("CS-SA01", "fastbact"), ("CS-SA01", "spades")],
            "runs follow dataset then harness tool order",
        )
        winner = runs[1]
        check(winner["assembly_sha256"] == "winner", "later metrics.json wins the merge_runs key")
        flat = harness.flatten_row(winner)
        check(list(flat) == list(harness.SUMMARY_COLUMNS), "flatten_row columns match SUMMARY_COLUMNS")
        with tsv_path.open(newline="", encoding="utf-8") as handle:
            tsv_rows = list(csv.reader(handle, delimiter="\t"))
        check(tsv_rows[0] == list(harness.SUMMARY_COLUMNS), "summary.tsv header matches schema")
        winner_cells = tsv_rows[2]
        check(len(winner_cells) == len(harness.SUMMARY_COLUMNS), "summary.tsv row width matches schema")
        expected_cells = ["" if value is None else str(value) for value in flat.values()]
        check(winner_cells == expected_cells, "TSV row matches flatten_row")
        check(winner_cells[list(harness.SUMMARY_COLUMNS).index("nga50")] == "", "null nga50 is an empty TSV cell")
        check(
            winner_cells[list(harness.SUMMARY_COLUMNS).index("command")] == "fastbact --example",
            "command list is flattened with the harness join",
        )

        dry_out = tmp / "dry-out"
        dry_code = main(["--raw-dir", str(raw), "--outdir", str(dry_out), "--dry-run"])
        check(dry_code == 0 and not dry_out.exists(), "dry-run validates without writing")

        empty = tmp / "empty"
        empty.mkdir()
        empty_out = tmp / "empty-out"
        empty_code = main(["--raw-dir", str(empty), "--outdir", str(empty_out)])
        empty_summary = json.loads((empty_out / "summary.json").read_text(encoding="utf-8"))
        empty_header = (empty_out / "summary.tsv").read_text(encoding="utf-8").splitlines()
        check(empty_code == 0 and empty_summary["runs"] == [], "empty raw dir writes zero runs")
        check(empty_header == ["\t".join(harness.SUMMARY_COLUMNS)], "empty TSV keeps the schema header")

        bad = tmp / "bad"
        _write_metrics(bad / "broken" / "metrics.json", {"dataset_id": "CS-SA01"})
        bad_code = main(["--raw-dir", str(bad), "--outdir", str(tmp / "bad-out")])
        check(bad_code == 2 and not (tmp / "bad-out").exists(), "incomplete metrics.json does not write a summary")

        broken = tmp / "broken-json"
        (broken / "run").mkdir(parents=True)
        (broken / "run" / "metrics.json").write_text("{", encoding="utf-8")
        broken_code = main(["--raw-dir", str(broken), "--outdir", str(tmp / "broken-out")])
        check(broken_code == 2, "invalid JSON is rejected")

    try:
        parse_args([])
    except SystemExit as exc:
        check(exc.code == 2, "raw-dir and outdir are required")
    else:
        check(False, "raw-dir and outdir are required")

    if failures:
        print(f"SELF_TEST_FAIL {len(failures)}")
        return 1
    print("SELF_TEST_PASS")
    return 0


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collate metrics.json files into summary.json and summary.tsv "
            "using the F26 competitor-benchmark schema."
        )
    )
    parser.add_argument("--raw-dir", type=Path, help="Directory searched recursively for metrics.json")
    parser.add_argument("--outdir", type=Path, help="Directory for summary.json and summary.tsv")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Load and validate metrics, print the planned summary, and write nothing",
    )
    parser.add_argument("--self-test", action="store_true", help="Run the built-in collation schema test")
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.self_test:
        return args
    missing = [name for name in ("raw_dir", "outdir") if getattr(args, name) is None]
    if missing:
        parser.error("required unless --self-test: " + ", ".join("--" + name.replace("_", "-") for name in missing))
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        return self_test()
    try:
        return collate(args.raw_dir, args.outdir, bool(args.dry_run))
    except CollationError as exc:
        print(f"[collate] {exc}", file=sys.stderr, flush=True)
        return 2


if __name__ == "__main__":
    sys.exit(main())
