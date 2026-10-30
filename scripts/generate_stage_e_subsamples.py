#!/usr/bin/env python3
"""Build the fastBact F26 Stage E depth-gradient FASTQs and manifest.

Reads a dataset manifest, calls scripts/subsample_fastq_pairs.py once per
dataset and target depth, and writes F26_STAGE_E_MANIFEST.tsv. The output
keeps every input column, so it stays compatible with
scripts/run_competitor_benchmark.py. fastq_r1 and fastq_r2 are set to the
subsampled files. Those columns override the harness read-discovery paths.

The on-disk F26_DATASET_MANIFEST.tsv columns are the schema that is written.
Dispatcher names map onto that schema as follows and are not renamed:

    organism -> species
    gram_status -> gram_stain
    gc_content -> gc_content_pct
    reference_fasta -> reference_fasta_path
    reference_sha256 -> reference_raw_sha256
    nominal_genome_size -> genome_size_bp
    platform -> instrument
    sra_run_accession -> sra_accession

dataset_id values look like CS-BS01_30x. A notes column records depth, seed,
source dataset, and the genome size used for the depth formula.

Genome size comes from genome_size_bp, nominal_genome_size, or genome_size
when that cell is a positive integer. Otherwise it is the base count of the
reference FASTA. Complete existing outputs are skipped unless --overwrite
is set. An output pair is complete only when both files exist and are
non-empty.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import math
import os
import random
import shlex
import subprocess
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SUBSAMPLE_SCRIPT = Path(__file__).resolve().parent / "subsample_fastq_pairs.py"
HARNESS_SCRIPT = Path(__file__).resolve().parent / "run_competitor_benchmark.py"

DEFAULT_INPUT_MANIFEST = REPO_ROOT / "F26_DATASET_MANIFEST.tsv"
DEFAULT_OUTPUT_MANIFEST = REPO_ROOT / "F26_STAGE_E_MANIFEST.tsv"
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "f26_cross_species" / "subsamples"
DEFAULT_DATA_ROOT = Path(
    "/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/04_ProkaDiff/data"
)
DEFAULT_DEPTHS = (30.0, 50.0, 100.0, 150.0)
DEFAULT_SEED = 42

# Kept in lockstep with run_competitor_benchmark.REQUIRED_MANIFEST_COLUMNS.
# The self-test imports the harness and fails if the tuples diverge.
HARNESS_REQUIRED_COLUMNS = (
    "dataset_id",
    "sra_accession",
    "bioproject",
    "reference_accession",
    "reference_fasta_path",
    "reference_raw_sha256",
)

SIZE_COLUMNS = ("genome_size_bp", "nominal_genome_size", "genome_size")
FASTA_COLUMNS = ("reference_fasta_path", "reference_fasta")
ACCESSION_COLUMNS = ("sra_accession", "sra_run_accession")
HISTORICAL_BIOPROJECT = "PRJNA1088182"


class StageEError(Exception):
    """A manifest, planning, or subsampling failure."""

    def __init__(self, message: str, exit_code: int = 1) -> None:
        super().__init__(message)
        self.exit_code = exit_code


def _positive_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and > 0")
    return value


def format_depth(depth: float) -> str:
    """Render 30.0 as 30 and 30.5 as 30.5, without scientific notation."""
    if not math.isfinite(depth) or depth <= 0.0:
        raise StageEError(f"depth must be finite and > 0, got {depth}")
    if float(depth).is_integer():
        return str(int(depth))
    return format(float(depth), "f").rstrip("0").rstrip(".")


def parse_genome_size(text: str) -> int | None:
    """Return a positive integer, or None when the cell is empty or not integral."""
    token = text.strip().replace(",", "")
    if not token:
        return None
    try:
        value = int(token)
    except ValueError:
        try:
            number = float(token)
        except ValueError:
            return None
        if not math.isfinite(number) or not number.is_integer():
            return None
        value = int(number)
    if value < 1:
        return None
    return value


def _field(row: Mapping[str, str], names: Sequence[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value.strip():
            return value.strip()
    return ""


def _first_present(fieldnames: Sequence[str], names: Sequence[str]) -> str | None:
    for name in names:
        if name in fieldnames:
            return name
    return None


def load_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Load a TSV manifest. Blank dataset_id rows are ignored."""
    if not path.is_file():
        raise StageEError(f"input manifest not found: {path}", exit_code=2)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise StageEError(f"manifest has no header: {path}", exit_code=2)
        fieldnames = [name for name in reader.fieldnames if name]
        missing = [name for name in HARNESS_REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise StageEError(
                f"manifest missing columns required by the benchmark harness {missing}: {path}",
                exit_code=2,
            )
        rows: list[dict[str, str]] = []
        seen: dict[str, int] = {}
        for row in reader:
            cleaned = {key: (value or "").strip() for key, value in row.items() if key}
            dataset_id = cleaned.get("dataset_id", "")
            if not dataset_id:
                continue
            seen[dataset_id] = seen.get(dataset_id, 0) + 1
            rows.append(cleaned)
    duplicates = sorted(dataset_id for dataset_id, count in seen.items() if count > 1)
    if duplicates:
        raise StageEError(
            f"duplicate dataset_id values in {path}: {', '.join(duplicates)}",
            exit_code=2,
        )
    if not rows:
        raise StageEError(f"manifest has no datasets: {path}", exit_code=2)
    return fieldnames, rows


def output_fieldnames(input_fieldnames: Sequence[str]) -> list[str]:
    """Input columns, then fastq paths, genome size, and notes when absent."""
    fields = list(input_fieldnames)
    for name in ("fastq_r1", "fastq_r2", "genome_size_bp", "notes"):
        if name not in fields:
            fields.append(name)
    return fields


def manifest_path(path: Path) -> str:
    """Repo-relative POSIX path when the logical path is inside the worktree.

    Symlinks are kept. ``data/f26_cross_species`` points at the shared data
    tree, and resolving it would record a path outside this worktree even
    though the harness opens ``data/...`` from the repository root.
    """
    logical = Path(os.path.abspath(path))
    try:
        return logical.relative_to(Path(os.path.abspath(REPO_ROOT))).as_posix()
    except ValueError:
        try:
            return str(path.resolve())
        except OSError:
            return str(logical)


def resolve_existing_file(raw: str, data_roots: Sequence[Path]) -> Path | None:
    """Find a manifest path against the cwd, the worktree, and known data roots."""
    if not raw:
        return None
    path = Path(raw)
    candidates: list[Path] = []
    if path.is_absolute():
        candidates.append(path)
    else:
        candidates.append(Path.cwd() / path)
        candidates.append(REPO_ROOT / path)
        stripped = Path(*path.parts[1:]) if path.parts and path.parts[0] == "data" else None
        for root in data_roots:
            if stripped is not None:
                candidates.append(root / stripped)
            candidates.append(root / path)
    seen: set[Path] = set()
    for candidate in candidates:
        try:
            key = candidate.resolve()
        except OSError:
            key = candidate
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            return candidate.resolve()
    return None


def genome_size_from_fasta(path: Path) -> int:
    """Count sequence characters in a FASTA file, summing every record."""
    opener = gzip_open_text if path.suffix == ".gz" or path.name.endswith(".gz") else None
    total = 0
    if opener is not None:
        handle_cm = opener(path)
    else:
        handle_cm = path.open("r", encoding="utf-8", newline="")
    with handle_cm as handle:
        for line in handle:
            if not line or line[0] == ">":
                continue
            total += sum(1 for char in line if char not in " \t\r\n")
    if total < 1:
        raise StageEError(f"reference FASTA contains no sequence bases: {path}")
    return total


def gzip_open_text(path: Path):
    import gzip

    return gzip.open(path, "rt", encoding="utf-8", newline="")


def default_data_roots() -> list[Path]:
    roots = [REPO_ROOT / "data", DEFAULT_DATA_ROOT]
    unique: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        try:
            key = root.resolve()
        except OSError:
            key = root
        if key in seen:
            continue
        seen.add(key)
        unique.append(root)
    return unique


def _read_name_variants(path: Path) -> list[Path]:
    text = str(path)
    if text.endswith(".fastq.gz"):
        return [path, Path(text[: -len(".gz")])]
    if text.endswith(".fq.gz"):
        return [path, Path(text[: -len(".gz")])]
    return [path]


def candidate_read_pairs(
    accession: str,
    bioproject: str,
    data_roots: Sequence[Path],
) -> list[tuple[Path, Path]]:
    """Same search order as run_competitor_benchmark.resolve_reads, plus plain FASTQ."""
    pairs: list[tuple[Path, Path]] = []
    for root in data_roots:
        historical = (
            root / "prjna1088182" / "fastq" / accession / f"{accession}_1.fastq.gz",
            root / "prjna1088182" / "fastq" / accession / f"{accession}_2.fastq.gz",
        )
        cross = (
            root / "f26_cross_species" / "fastq" / f"{accession}_1.fastq.gz",
            root / "f26_cross_species" / "fastq" / f"{accession}_2.fastq.gz",
        )
        ordered = (historical, cross) if bioproject == HISTORICAL_BIOPROJECT else (cross, historical)
        for left, right in ordered:
            for r1, r2 in zip(_read_name_variants(left), _read_name_variants(right)):
                pairs.append((r1, r2))
    return pairs


def resolve_reads(row: Mapping[str, str], data_roots: Sequence[Path]) -> tuple[Path, Path]:
    """Use explicit fastq columns when both are set. Otherwise discover by accession."""
    raw_r1 = row.get("fastq_r1", "").strip()
    raw_r2 = row.get("fastq_r2", "").strip()
    dataset_id = row.get("dataset_id", "")
    if raw_r1 or raw_r2:
        if not raw_r1 or not raw_r2:
            raise StageEError(
                f"{dataset_id}: fastq_r1 and fastq_r2 must both be set when either is set"
            )
        r1 = resolve_existing_file(raw_r1, data_roots)
        r2 = resolve_existing_file(raw_r2, data_roots)
        if r1 is None or r2 is None:
            missing = raw_r1 if r1 is None else raw_r2
            raise StageEError(f"{dataset_id}: FASTQ not found: {missing}")
        return r1, r2

    accession = _field(row, ACCESSION_COLUMNS)
    if not accession:
        raise StageEError(f"{dataset_id}: no fastq_r1/fastq_r2 and no sra_accession")
    bioproject = row.get("bioproject", "").strip()
    tried = candidate_read_pairs(accession, bioproject, data_roots)
    for r1, r2 in tried:
        if r1.is_file() and r2.is_file():
            return r1.resolve(), r2.resolve()
    preview = ", ".join(str(r1) for r1, _r2 in tried[:4])
    raise StageEError(f"{dataset_id}: FASTQ pair not found for {accession}. Tried: {preview}")


def resolve_genome_size(
    row: Mapping[str, str],
    fieldnames: Sequence[str],
    data_roots: Sequence[Path],
) -> tuple[int, str | None]:
    """Return (genome size, column to fill when the tabulated cell was empty)."""
    for name in SIZE_COLUMNS:
        if name not in row:
            continue
        parsed = parse_genome_size(row.get(name, ""))
        if parsed is not None:
            return parsed, None
    present = _first_present(fieldnames, SIZE_COLUMNS)
    raw_fasta = _field(row, FASTA_COLUMNS)
    dataset_id = row.get("dataset_id", "")
    if not raw_fasta:
        raise StageEError(
            f"{dataset_id}: genome size is not tabulated and reference_fasta_path is empty"
        )
    fasta = resolve_existing_file(raw_fasta, data_roots)
    if fasta is None:
        raise StageEError(
            f"{dataset_id}: genome size is not tabulated and reference FASTA was not found: {raw_fasta}"
        )
    return genome_size_from_fasta(fasta), present or "genome_size_bp"


def _check_dataset_token(dataset_id: str) -> None:
    if (
        not dataset_id
        or dataset_id in {".", ".."}
        or "/" in dataset_id
        or "\\" in dataset_id
        or "\x00" in dataset_id
    ):
        raise StageEError(f"dataset_id cannot be used as a path component: {dataset_id!r}", exit_code=2)


def _same_path(left: Path, right: Path) -> bool:
    try:
        return left.resolve() == right.resolve()
    except OSError:
        return os.path.abspath(left) == os.path.abspath(right)


def subsample_command(
    r1: Path,
    r2: Path,
    out1: Path,
    out2: Path,
    depth: float,
    genome_size: int,
    seed: int,
) -> list[str]:
    return [
        sys.executable,
        str(SUBSAMPLE_SCRIPT),
        "--r1",
        str(r1),
        "--r2",
        str(r2),
        "--out1",
        str(out1),
        "--out2",
        str(out2),
        "--target-depth",
        format_depth(depth),
        "--genome-size",
        str(genome_size),
        "--seed",
        str(seed),
    ]


def outputs_complete(out1: Path, out2: Path) -> bool:
    try:
        return out1.is_file() and out2.is_file() and out1.stat().st_size > 0 and out2.stat().st_size > 0
    except OSError:
        return False


class PlannedJob:
    def __init__(
        self,
        source_id: str,
        dataset_id: str,
        depth: float,
        genome_size: int,
        r1: Path,
        r2: Path,
        out1: Path,
        out2: Path,
        row: dict[str, str],
        command: list[str],
    ) -> None:
        self.source_id = source_id
        self.dataset_id = dataset_id
        self.depth = depth
        self.genome_size = genome_size
        self.r1 = r1
        self.r2 = r2
        self.out1 = out1
        self.out2 = out2
        self.row = row
        self.command = command

    @property
    def depth_label(self) -> str:
        return format_depth(self.depth)


def stage_note(source_id: str, depth_label: str, seed: int, genome_size: int, previous: str) -> str:
    note = (
        f"stage_e depth={depth_label}x seed={seed} source={source_id} "
        f"genome_size_bp={genome_size}"
    )
    if previous:
        return previous + " | " + note
    return note


def plan_jobs(
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[str, str]],
    depths: Sequence[float],
    out_dir: Path,
    seed: int,
    data_roots: Sequence[Path],
) -> list[PlannedJob]:
    if not SUBSAMPLE_SCRIPT.is_file():
        raise StageEError(f"subsampler not found: {SUBSAMPLE_SCRIPT}", exit_code=2)
    labels = [format_depth(depth) for depth in depths]
    if len(labels) != len(set(labels)):
        raise StageEError(f"duplicate depths after formatting: {', '.join(labels)}", exit_code=2)

    out_fields = output_fieldnames(fieldnames)
    jobs: list[PlannedJob] = []
    for source in rows:
        source_id = source["dataset_id"]
        _check_dataset_token(source_id)
        r1, r2 = resolve_reads(source, data_roots)
        genome_size, size_column = resolve_genome_size(source, fieldnames, data_roots)
        for depth, label in zip(depths, labels):
            dataset_id = f"{source_id}_{label}x"
            _check_dataset_token(dataset_id)
            sample_dir = out_dir / dataset_id
            out1 = sample_dir / f"{dataset_id}_R1.fastq.gz"
            out2 = sample_dir / f"{dataset_id}_R2.fastq.gz"
            for source_path, dest in ((r1, out1), (r1, out2), (r2, out1), (r2, out2)):
                if _same_path(source_path, dest):
                    raise StageEError(
                        f"{dataset_id}: refusing to write subsample output over input {source_path}"
                    )
            row = {name: source.get(name, "") for name in out_fields}
            row["dataset_id"] = dataset_id
            row["fastq_r1"] = manifest_path(out1)
            row["fastq_r2"] = manifest_path(out2)
            if size_column is not None:
                row[size_column] = str(genome_size)
            elif "genome_size_bp" in row and not row["genome_size_bp"]:
                row["genome_size_bp"] = str(genome_size)
            row["notes"] = stage_note(source_id, label, seed, genome_size, source.get("notes", ""))
            command = subsample_command(r1, r2, out1, out2, depth, genome_size, seed)
            jobs.append(
                PlannedJob(
                    source_id=source_id,
                    dataset_id=dataset_id,
                    depth=depth,
                    genome_size=genome_size,
                    r1=r1,
                    r2=r2,
                    out1=out1,
                    out2=out2,
                    row=row,
                    command=command,
                )
            )
    return jobs


def render_manifest(fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(fieldnames),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") for name in fieldnames})
    return buffer.getvalue()


def write_manifest(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = render_manifest(fieldnames, rows)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(payload, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def execute(args: argparse.Namespace) -> int:
    fieldnames, rows = load_manifest(Path(args.input_manifest))
    if args.dataset:
        selected = [row for row in rows if row["dataset_id"] == args.dataset]
        if not selected:
            known = ", ".join(row["dataset_id"] for row in rows)
            raise StageEError(f"unknown dataset {args.dataset!r}. Known: {known}", exit_code=2)
        rows = selected

    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = (Path.cwd() / out_dir).resolve()
    jobs = plan_jobs(
        fieldnames,
        rows,
        list(args.depths),
        out_dir,
        int(args.seed),
        default_data_roots(),
    )
    out_fields = output_fieldnames(fieldnames)
    overwrite = bool(args.overwrite)

    if args.dry_run:
        for job in jobs:
            action = "skip" if outputs_complete(job.out1, job.out2) and not overwrite else "run"
            print(f"[stage-e] dry-run {action} {shlex.join(job.command)}", file=sys.stderr, flush=True)
        sys.stdout.write(render_manifest(out_fields, [job.row for job in jobs]))
        return 0

    generated = 0
    skipped = 0
    for job in jobs:
        if outputs_complete(job.out1, job.out2) and not overwrite:
            skipped += 1
            print(
                f"[stage-e] skip {job.dataset_id}: outputs exist and are non-empty",
                file=sys.stderr,
                flush=True,
            )
            continue
        job.out1.parent.mkdir(parents=True, exist_ok=True)
        print(f"[stage-e] run {shlex.join(job.command)}", file=sys.stderr, flush=True)
        result = subprocess.run(job.command, check=False)
        if result.returncode != 0:
            raise StageEError(
                f"subsampler exited {result.returncode} for {job.dataset_id}"
            )
        if not outputs_complete(job.out1, job.out2):
            raise StageEError(f"subsampler did not create non-empty outputs for {job.dataset_id}")
        generated += 1

    manifest = Path(args.output_manifest)
    write_manifest(manifest, out_fields, [job.row for job in jobs])
    print(
        f"[stage-e] done rows={len(jobs)} generated={generated} skipped={skipped} "
        f"manifest={manifest}",
        file=sys.stderr,
        flush=True,
    )
    return 0


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate F26 Stage E depth-gradient FASTQ pairs and a benchmark manifest. "
            "Each dataset is subsampled with scripts/subsample_fastq_pairs.py."
        ),
        epilog=(
            "examples:\n"
            "  python3 scripts/generate_stage_e_subsamples.py --dataset CS-BS01\n"
            "  python3 scripts/generate_stage_e_subsamples.py --dry-run\n"
            "  python3 scripts/generate_stage_e_subsamples.py --self-test\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST,
        help=f"Source dataset TSV (default: {DEFAULT_INPUT_MANIFEST.name})",
    )
    parser.add_argument(
        "--depths",
        type=_positive_float,
        nargs="+",
        default=list(DEFAULT_DEPTHS),
        help="Target depths (default: 30 50 100 150)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Directory for subsampled FASTQ pairs",
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=DEFAULT_OUTPUT_MANIFEST,
        help=f"Stage E manifest to write (default: {DEFAULT_OUTPUT_MANIFEST.name})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Sampling seed passed to the subsampler (default: {DEFAULT_SEED})",
    )
    parser.add_argument(
        "--dataset",
        help="Generate rows for this dataset_id only",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print subsampler commands and the planned manifest without writing FASTQs",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate outputs even when both FASTQs already exist and are non-empty",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip complete outputs (this is already the default; mutually exclusive with --overwrite)",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run the built-in manifest, execution, and determinism tests",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.overwrite and args.skip_existing:
        parser.error("--overwrite and --skip-existing are mutually exclusive")
    return args


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_fastq(path: Path, n_pairs: int, mate: int, read_len: int = 50) -> None:
    lines: list[str] = []
    pad = "A" if mate == 1 else "C"
    for index in range(n_pairs):
        token = f"{index:08d}"
        sequence = pad * (read_len - len(token)) + token
        lines.append(f"@pair{index:06d}/{mate}\n{sequence}\n+\n{'I' * read_len}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="ascii")


def _write_manifest(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    write_manifest(path, fieldnames, rows)


def _load_harness():
    spec = importlib.util.spec_from_file_location("f26_stage_e_harness_check", HARNESS_SCRIPT)
    if spec is None or spec.loader is None:
        raise AssertionError(f"cannot load harness: {HARNESS_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _capture_main(argv: Sequence[str]) -> tuple[int, str, str]:
    stdout = io.StringIO()
    stderr = io.StringIO()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = main(list(argv))
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
    return int(code), stdout.getvalue(), stderr.getvalue()


def self_test() -> int:
    """Execution, schema, skip/overwrite, determinism, and dry-run checks."""
    failures: list[str] = []

    def check(condition: bool, message: str) -> None:
        if not condition:
            failures.append(message)
            print(f"FAIL {message}", file=sys.stderr)

    linked = manifest_path(
        REPO_ROOT / "data" / "f26_cross_species" / "subsamples" / "CS-BS01_30x" / "CS-BS01_30x_R1.fastq.gz"
    )
    check(
        linked == "data/f26_cross_species/subsamples/CS-BS01_30x/CS-BS01_30x_R1.fastq.gz",
        f"output paths stay repo-relative across the data symlink ({linked})",
    )

    real_manifest = REPO_ROOT / "F26_DATASET_MANIFEST.tsv"
    real_before = real_manifest.read_bytes() if real_manifest.is_file() else None
    default_manifest = DEFAULT_OUTPUT_MANIFEST
    default_existed = default_manifest.exists()
    default_before = default_manifest.read_bytes() if default_existed else None

    try:
        harness = _load_harness()
        check(
            tuple(harness.REQUIRED_MANIFEST_COLUMNS) == HARNESS_REQUIRED_COLUMNS,
            "harness required columns match this script",
        )
    except Exception as exc:  # noqa: BLE001 - self-test must report loader failures
        check(False, f"harness import failed: {exc}")
        harness = None

    with tempfile.TemporaryDirectory(prefix="f26-stage-e-self-test-") as tmp_name:
        tmp = Path(tmp_name)
        n_pairs = 16
        read_len = 50
        genome_size = 40
        r1 = tmp / "reads" / "SYN_1.fastq"
        r2 = tmp / "reads" / "SYN_2.fastq"
        _write_fastq(r1, n_pairs, 1, read_len)
        _write_fastq(r2, n_pairs, 2, read_len)
        fasta = tmp / "ref.fasta"
        fasta.write_text(">chr1\nACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n>chr2\nAAAA\n", encoding="ascii")
        check(genome_size_from_fasta(fasta) == 44, "multi-record FASTA base count includes every record")
        short_fasta = tmp / "short.fasta"
        short_fasta.write_text(">chr\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="ascii")
        check(genome_size_from_fasta(short_fasta) == genome_size, "single-record FASTA base count")

        real_header = real_manifest.read_text(encoding="utf-8").splitlines()[0].split("\t")
        check("dataset_id" in real_header and "genome_size_bp" in real_header, "real manifest header is readable")
        fieldnames = list(real_header)
        if "fastq_r1" not in fieldnames:
            fieldnames.extend(["fastq_r1", "fastq_r2"])

        def source_row(dataset_id: str, accession: str, size: str, fasta_path: str) -> dict[str, str]:
            row = {name: "" for name in fieldnames}
            row.update(
                {
                    "dataset_id": dataset_id,
                    "species": "Syntheticomonas test",
                    "strain": "ST-1",
                    "gram_stain": "Negative",
                    "gc_content_pct": "50.00",
                    "genome_size_bp": size,
                    "sra_accession": accession,
                    "bioproject": "PRJNA0",
                    "biosample": "SAMN0",
                    "library_layout": "PAIRED (2x50)",
                    "instrument": "synthetic",
                    "total_spots": str(n_pairs),
                    "reference_accession": "SYN_REF",
                    "reference_fasta_path": fasta_path,
                    "reference_raw_sha256": "a" * 64,
                    "reference_normalized_seq_sha256": "b" * 64,
                    "role": "F26 Stage E self-test",
                    "fastq_r1": str(r1),
                    "fastq_r2": str(r2),
                }
            )
            return {name: row.get(name, "") for name in fieldnames}

        input_rows = [
            source_row("SYN-A", "SYNACC_A", str(genome_size), str(short_fasta)),
            source_row("SYN-B", "SYNACC_B", str(genome_size), "/does/not/exist.fasta"),
        ]
        input_manifest = tmp / "mini_manifest.tsv"
        _write_manifest(input_manifest, fieldnames, input_rows)
        input_hash = _sha256_file(input_manifest)

        out_dir = tmp / "subsamples"
        output_manifest = tmp / "F26_STAGE_E_MANIFEST.tsv"
        code, stdout, stderr = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "50",
                "--out-dir",
                str(out_dir),
                "--output-manifest",
                str(output_manifest),
                "--seed",
                "42",
                "--dataset",
                "SYN-A",
            ]
        )
        check(code == 0, f"filtered execution exits 0 (got {code}; stderr={stderr!r})")
        check(stdout == "", "execution does not write the manifest to stdout")
        check(_sha256_file(input_manifest) == input_hash, "input manifest bytes are unchanged")
        check(output_manifest.is_file(), "output manifest was written")

        with output_manifest.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            out_fields = list(reader.fieldnames or [])
            out_rows = list(reader)
        check(out_fields[: len(fieldnames)] == fieldnames, "output columns start with the input schema")
        for name in ("fastq_r1", "fastq_r2", "notes"):
            check(name in out_fields, f"output manifest has {name}")
        check([row["dataset_id"] for row in out_rows] == ["SYN-A_30x", "SYN-A_50x"], "dataset ids encode depth")
        check(all(row["species"] == "Syntheticomonas test" for row in out_rows), "organism column is preserved")
        check(
            all(row["reference_raw_sha256"] == "a" * 64 for row in out_rows),
            "reference hash is copied unchanged",
        )
        check(all(row["genome_size_bp"] == "40" for row in out_rows), "tabulated genome size is kept")
        check(all("source=SYN-A" in row["notes"] and "seed=42" in row["notes"] for row in out_rows), "notes record provenance")
        check("SYN-B" not in output_manifest.read_text(encoding="utf-8"), "--dataset restricts planning")

        if harness is not None:
            loaded = harness.load_manifest(output_manifest)
            check(len(loaded) == 2, "benchmark harness accepts the Stage E manifest")
            check(loaded[0]["dataset_id"] == "SYN-A_30x", "harness reads the depth dataset id")
            check(loaded[0]["fastq_r1"].endswith("SYN-A_30x_R1.fastq.gz"), "harness sees the subsampled R1 path")

        out_30_r1 = out_dir / "SYN-A_30x" / "SYN-A_30x_R1.fastq.gz"
        out_30_r2 = out_dir / "SYN-A_30x" / "SYN-A_30x_R2.fastq.gz"
        out_50_r1 = out_dir / "SYN-A_50x" / "SYN-A_50x_R1.fastq.gz"
        check(outputs_complete(out_30_r1, out_30_r2), "30x outputs exist and are non-empty")
        check(out_50_r1.is_file() and out_50_r1.stat().st_size > 0, "50x R1 exists and is non-empty")
        hash_30 = _sha256_file(out_30_r1)
        hash_30_r2 = _sha256_file(out_30_r2)
        check(hash_30 != hash_30_r2, "R1 and R2 subsamples are not copies")
        check(out_30_r1.stat().st_size < r1.stat().st_size, "30x output is smaller than the input FASTQ")
        expected_k = int(round((30 * genome_size) / (2 * read_len)))
        check(expected_k == 12 and expected_k < n_pairs, "self-test depth selects a proper subset")

        code_skip, _stdout_skip, stderr_skip = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(out_dir),
                "--output-manifest",
                str(tmp / "manifest_skip.tsv"),
                "--dataset",
                "SYN-A",
                "--skip-existing",
            ]
        )
        check(code_skip == 0, f"skip-existing exits 0 (got {code_skip})")
        check("skip SYN-A_30x" in stderr_skip, "skip-existing reports the completed pair")
        check(_sha256_file(out_30_r1) == hash_30, "skip-existing does not rewrite a complete R1")

        other_dir = tmp / "subsamples_repeat"
        code_repeat, _stdout_repeat, _stderr_repeat = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(other_dir),
                "--output-manifest",
                str(tmp / "manifest_repeat.tsv"),
                "--seed",
                "42",
                "--dataset",
                "SYN-A",
            ]
        )
        repeated = other_dir / "SYN-A_30x" / "SYN-A_30x_R1.fastq.gz"
        check(code_repeat == 0 and repeated.is_file(), "repeat run with the same seed exits 0")
        check(_sha256_file(repeated) == hash_30, "same seed reproduces the R1 SHA-256")

        other_seed_dir = tmp / "subsamples_seed99"
        code_seed, _stdout_seed, _stderr_seed = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(other_seed_dir),
                "--output-manifest",
                str(tmp / "manifest_seed99.tsv"),
                "--seed",
                "99",
                "--dataset",
                "SYN-A",
            ]
        )
        seeded = other_seed_dir / "SYN-A_30x" / "SYN-A_30x_R1.fastq.gz"
        check(code_seed == 0 and seeded.is_file(), "alternate seed exits 0")
        check(_sha256_file(seeded) != hash_30, "a different seed changes the subsample")
        oracle_42 = sorted(random.Random(42).sample(range(n_pairs), expected_k))
        oracle_99 = sorted(random.Random(99).sample(range(n_pairs), expected_k))
        check(oracle_42 != oracle_99, "self-test seeds select different pair indices")

        out_30_r1.write_bytes(b"STALE")
        out_30_r2.write_bytes(b"STALE")
        code_keep, _stdout_keep, _stderr_keep = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(out_dir),
                "--output-manifest",
                str(tmp / "manifest_keep.tsv"),
                "--dataset",
                "SYN-A",
            ]
        )
        check(code_keep == 0 and out_30_r1.read_bytes() == b"STALE", "non-empty existing outputs are kept")
        code_overwrite, _stdout_overwrite, stderr_overwrite = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(out_dir),
                "--output-manifest",
                str(tmp / "manifest_overwrite.tsv"),
                "--dataset",
                "SYN-A",
                "--overwrite",
            ]
        )
        check(code_overwrite == 0, f"overwrite exits 0 (got {code_overwrite}; stderr={stderr_overwrite!r})")
        check(out_30_r1.read_bytes()[:2] == b"\x1f\x8b", "overwrite replaces stale bytes with gzip")
        check(_sha256_file(out_30_r1) == hash_30, "overwrite with the same seed restores the original digest")

        empty_r1 = out_dir / "SYN-A_30x" / "SYN-A_30x_R1.fastq.gz"
        empty_r1.write_bytes(b"")
        code_empty, _stdout_empty, stderr_empty = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(out_dir),
                "--output-manifest",
                str(tmp / "manifest_empty.tsv"),
                "--dataset",
                "SYN-A",
            ]
        )
        check(code_empty == 0 and "run " in stderr_empty, "an empty output is regenerated")
        check(empty_r1.stat().st_size > 0, "regenerated output is non-empty")

        dry_dir = tmp / "dry"
        code_dry, stdout_dry, stderr_dry = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(dry_dir),
                "--output-manifest",
                str(tmp / "should_not_be_written.tsv"),
                "--dataset",
                "SYN-A",
                "--dry-run",
            ]
        )
        check(code_dry == 0, f"dry-run exits 0 (got {code_dry})")
        check(not dry_dir.exists(), "dry-run does not create the output directory")
        check(not (tmp / "should_not_be_written.tsv").exists(), "dry-run does not write the manifest file")
        check("--target-depth 30" in stderr_dry and "--genome-size 40" in stderr_dry and "--seed 42" in stderr_dry, "dry-run prints the subsampler command")
        check("subsample_fastq_pairs.py" in stderr_dry, "dry-run names the subsampler")
        dry_reader = csv.DictReader(io.StringIO(stdout_dry), delimiter="\t")
        dry_rows = list(dry_reader)
        check(
            [row["dataset_id"] for row in dry_rows] == ["SYN-A_30x"],
            "dry-run stdout is the planned manifest",
        )

        missing_size_fields = [name for name in fieldnames if name != "genome_size_bp"]
        missing_size_row = source_row("SYN-C", "SYNACC_C", "", str(short_fasta))
        missing_size_row = {name: missing_size_row.get(name, "") for name in missing_size_fields}
        missing_manifest = tmp / "missing_size.tsv"
        _write_manifest(missing_manifest, missing_size_fields, [missing_size_row])
        code_size, _stdout_size, stderr_size = _capture_main(
            [
                "--input-manifest",
                str(missing_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(tmp / "size_out"),
                "--output-manifest",
                str(tmp / "size_manifest.tsv"),
                "--dataset",
                "SYN-C",
                "--dry-run",
            ]
        )
        check(code_size == 0, f"missing tabulated size still plans (got {code_size}; stderr={stderr_size!r})")
        check("--genome-size 40" in stderr_size, "genome size is calculated from the reference FASTA")

        discover_root = tmp / "data"
        accession = "SYNACC_D"
        discovered_r1 = discover_root / "f26_cross_species" / "fastq" / f"{accession}_1.fastq"
        discovered_r2 = discover_root / "f26_cross_species" / "fastq" / f"{accession}_2.fastq"
        _write_fastq(discovered_r1, 2, 1, 20)
        _write_fastq(discovered_r2, 2, 2, 20)
        found_r1, found_r2 = resolve_reads(
            {"dataset_id": "SYN-D", "sra_accession": accession, "bioproject": "PRJNA1"},
            [discover_root],
        )
        check(found_r1 == discovered_r1.resolve() and found_r2 == discovered_r2.resolve(), "accession discovery finds the cross-species pair")

        clash_r1 = tmp / "clash" / "SYN-E_30x" / "SYN-E_30x_R1.fastq.gz"
        clash_r1.parent.mkdir(parents=True)
        clash_r1.write_bytes(b"INPUT")
        clash_row = source_row("SYN-E", "SYNACC_E", "40", str(short_fasta))
        clash_row["fastq_r1"] = str(clash_r1)
        clash_row["fastq_r2"] = str(r2)
        clash_manifest = tmp / "clash.tsv"
        _write_manifest(clash_manifest, fieldnames, [clash_row])
        code_clash, _stdout_clash, stderr_clash = _capture_main(
            [
                "--input-manifest",
                str(clash_manifest),
                "--depths",
                "30",
                "--out-dir",
                str(tmp / "clash"),
                "--output-manifest",
                str(tmp / "clash_out.tsv"),
                "--dataset",
                "SYN-E",
            ]
        )
        check(code_clash == 1 and clash_r1.read_bytes() == b"INPUT", "an output path that aliases an input is refused")
        check("refusing" in stderr_clash, "the refusal is reported")

        code_unknown, _stdout_unknown, stderr_unknown = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--dataset",
                "MISSING",
                "--out-dir",
                str(tmp / "unused"),
                "--output-manifest",
                str(tmp / "unused.tsv"),
                "--dry-run",
            ]
        )
        check(code_unknown == 2 and "unknown dataset" in stderr_unknown, "unknown --dataset exits 2")

        code_flags, _stdout_flags, stderr_flags = _capture_main(
            ["--overwrite", "--skip-existing", "--dry-run"]
        )
        check(code_flags == 2 and "mutually exclusive" in stderr_flags, "overwrite and skip-existing conflict exits 2")

        code_missing_fasta, _stdout_missing, stderr_missing = _capture_main(
            [
                "--input-manifest",
                str(input_manifest),
                "--dataset",
                "SYN-B",
                "--depths",
                "30",
                "--out-dir",
                str(tmp / "missing_fasta_out"),
                "--output-manifest",
                str(tmp / "missing_fasta.tsv"),
            ]
        )
        check(
            code_missing_fasta == 0 and (tmp / "missing_fasta_out" / "SYN-B_30x" / "SYN-B_30x_R1.fastq.gz").is_file(),
            "a tabulated genome size does not require the reference FASTA to exist",
        )
        check("error:" not in stderr_missing, "tabulated-size run does not report an error")

    if real_before is not None:
        check(real_manifest.read_bytes() == real_before, "F26_DATASET_MANIFEST.tsv is not modified")
    if default_existed:
        check(default_manifest.read_bytes() == default_before, "default Stage E manifest is not modified")
    else:
        check(not default_manifest.exists(), "self-test does not create the default Stage E manifest")

    if failures:
        print(f"SELF_TEST_FAIL {len(failures)}", file=sys.stderr)
        return 1
    print("SELF_TEST_PASS")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args = parse_args(sys.argv[1:] if argv is None else argv)
    except SystemExit:
        raise
    if args.self_test:
        return self_test()
    try:
        return execute(args)
    except StageEError as exc:
        print(f"error: {exc}", file=sys.stderr, flush=True)
        return exc.exit_code
    except OSError as exc:
        print(f"error: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
