#!/usr/bin/env python3
"""Deterministic paired-end FASTQ subsampler for fastBact Phase F26 Stage E.

Draws a reproducible subset of read pairs for coverage-sensitivity benchmarks
(for example 30x, 50x, 100x, and 150x) without loading either FASTQ into memory.

Algorithm
    Pass 1 streams R1 and R2 in lockstep, counts synchronized pairs N, checks
    mate headers, and reads the first R1 length when --read-len is omitted.
    Pass 2 streams both files again and writes the selected pairs. Indices are
    ``sorted(random.Random(seed).sample(range(N), k))``, so mates stay together
    and the output order is the original file order.

Depth formula
    k = round((target_depth * genome_size) / (2 * read_len))
    Python round (half to even) is used. If k is greater than N, every pair is
    kept and the run still exits 0.

Gzip outputs use a fixed mtime of 0 and omit the original filename so the
file SHA-256 depends on the reads, the seed, and the compression level, not
on the clock or the output path.

Standard library only. Records are 4-line FASTQ (Illumina / SRA layout).
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import io
import random
import sys
import tempfile
import zlib
from pathlib import Path
from typing import BinaryIO, Iterator, Sequence

GZIP_COMPRESSLEVEL = 6
GZIP_MTIME = 0
DEFAULT_SEED = 42


class SubsampleError(Exception):
    """A user-facing subsampling failure."""


class FastqFormatError(SubsampleError):
    """A FASTQ record is truncated or malformed."""


class PairSyncError(SubsampleError):
    """R1 and R2 are not a synchronized pair stream."""


def pairs_for_depth(depth: float, genome_size: int, read_len: int) -> int:
    """Pair count for a target depth. k = round((depth * genome) / (2 * read_len))."""
    if read_len < 1:
        raise SubsampleError(f"read length must be >= 1, got {read_len}")
    if depth <= 0 or genome_size < 1:
        raise SubsampleError("target depth and genome size must be positive")
    return int(round((depth * genome_size) / (2 * read_len)))


def pairs_for_fraction(fraction: float, n_pairs: int) -> int:
    """Pair count for a fraction of N. Uses the same half-to-even round."""
    if not 0.0 < fraction <= 1.0:
        raise SubsampleError(f"fraction must satisfy 0 < f <= 1, got {fraction}")
    if n_pairs < 0:
        raise SubsampleError(f"pair count must be >= 0, got {n_pairs}")
    return int(round(fraction * n_pairs))


def select_pair_indices(n_pairs: int, k: int, seed: int) -> list[int] | None:
    """Sorted 0-based indices from Random(seed).sample(range(N), k).

    Returns None when k == N. Sorting a sample of the entire population is
    exactly range(N) for every seed, so the keep-all path does not allocate
    a permutation of the file.
    """
    if k < 0 or n_pairs < 0 or k > n_pairs:
        raise SubsampleError(f"cannot select {k} pairs from {n_pairs}")
    if k == n_pairs:
        return None
    rng = random.Random(seed)
    return sorted(rng.sample(range(n_pairs), k))


def _payload(line: bytes) -> bytes:
    if line.endswith(b"\r\n"):
        return line[:-2]
    if line.endswith((b"\n", b"\r")):
        return line[:-1]
    return line


def _illumina_space_mate(comment: bytes) -> int | None:
    """Mate number from an Illumina comment tag such as ``1:N:0:ACGT``."""
    if len(comment) >= 2 and comment[:1] in (b"1", b"2") and comment[1:2] == b":":
        return 1 if comment[:1] == b"1" else 2
    return None


def split_header(header_line: bytes) -> tuple[bytes, int | None]:
    """Return (canonical read id, mate number or None).

    ``/1`` and ``/2`` suffixes on the first token, and Illumina space tags
    (``1:...`` / ``2:...``), are recognized as mate labels. The canonical id
    is what remains after the suffix is removed, so ``@id/1`` matches ``@id/2``
    and identical SRA headers match each other.
    """
    line = header_line.rstrip(b"\r\n")
    if not line.startswith(b"@"):
        raise FastqFormatError("FASTQ header does not start with '@'")
    body = line[1:]
    if not body:
        raise FastqFormatError("FASTQ header has an empty read id")
    parts = body.split(None, 1)
    token = parts[0]
    comment = parts[1] if len(parts) > 1 else b""
    mate: int | None = None
    if token.endswith(b"/1") or token.endswith(b"/2"):
        mate = 1 if token.endswith(b"/1") else 2
        token = token[:-2]
    if not token:
        raise FastqFormatError("empty read id after removing the /1 or /2 suffix")
    tag = _illumina_space_mate(comment)
    if tag is not None:
        if mate is not None and mate != tag:
            shown = line.decode("utf-8", "replace")
            raise FastqFormatError(
                f"mate suffix /{mate} conflicts with Illumina tag {tag} in {shown}"
            )
        mate = tag
    return token, mate


def assert_pair_headers(header_r1: bytes, header_r2: bytes, pair_number: int) -> None:
    """Require R1 and R2 to name the same spot, with mate order 1 then 2 when both are labeled."""
    try:
        id1, mate1 = split_header(header_r1)
        id2, mate2 = split_header(header_r2)
    except FastqFormatError as exc:
        raise FastqFormatError(f"pair {pair_number}: {exc}") from exc
    if id1 != id2:
        left = id1.decode("utf-8", "replace")
        right = id2.decode("utf-8", "replace")
        raise PairSyncError(
            f"pair {pair_number}: R1/R2 read ids differ ({left} vs {right})"
        )
    if mate1 is None or mate2 is None or mate1 == mate2:
        return
    if mate1 == 1 and mate2 == 2:
        return
    raise PairSyncError(
        f"pair {pair_number}: mate tags are R1={mate1} R2={mate2}; expected 1 then 2"
    )


def open_fastq_input(path: Path) -> BinaryIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_record(
    handle: BinaryIO,
    path: Path,
    pair_number: int,
) -> tuple[bytes, bytes, bytes, bytes] | None:
    """Read one 4-line FASTQ record, or None at a clean EOF."""
    header = handle.readline()
    if header == b"":
        return None
    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()
    if seq == b"" or plus == b"" or qual == b"":
        raise FastqFormatError(f"{path}: truncated FASTQ record at pair {pair_number}")
    if not header.startswith(b"@"):
        raise FastqFormatError(
            f"{path}: pair {pair_number} header does not start with '@'"
        )
    if not plus.startswith(b"+"):
        raise FastqFormatError(
            f"{path}: pair {pair_number} separator does not start with '+'"
        )
    seq_len = len(_payload(seq))
    qual_len = len(_payload(qual))
    if seq_len != qual_len:
        raise FastqFormatError(
            f"{path}: pair {pair_number} sequence length {seq_len} "
            f"differs from quality length {qual_len}"
        )
    if seq_len == 0:
        raise FastqFormatError(f"{path}: pair {pair_number} has an empty sequence")
    return header, seq, plus, qual


def iter_pair_records(
    r1_path: Path,
    r2_path: Path,
) -> Iterator[tuple[int, tuple[bytes, bytes, bytes, bytes], tuple[bytes, bytes, bytes, bytes]]]:
    """Yield (1-based pair number, R1 record, R2 record) until both files end together."""
    with open_fastq_input(r1_path) as r1_handle, open_fastq_input(r2_path) as r2_handle:
        pair_number = 0
        while True:
            pair_number += 1
            left = read_record(r1_handle, r1_path, pair_number)
            right = read_record(r2_handle, r2_path, pair_number)
            if left is None and right is None:
                return
            if left is None or right is None:
                ended = "R1" if left is None else "R2"
                raise PairSyncError(
                    f"R1 and R2 pair counts differ at pair {pair_number} ({ended} ended first)"
                )
            assert_pair_headers(left[0], right[0], pair_number)
            yield pair_number, left, right


def scan_pairs(r1_path: Path, r2_path: Path) -> tuple[int, int, bool]:
    """Pass 1: count pairs, take the first R1 length, and validate every header.

    Returns (n_pairs, first_r1_read_length, variable_r1_length).
    """
    n_pairs = 0
    first_len: int | None = None
    variable = False
    for _pair_number, left, _right in iter_pair_records(r1_path, r2_path):
        seq_len = len(_payload(left[1]))
        if first_len is None:
            first_len = seq_len
        elif seq_len != first_len:
            variable = True
        n_pairs += 1
    if n_pairs == 0 or first_len is None:
        raise FastqFormatError("input contains no read pairs")
    return n_pairs, first_len, variable


class _HashingSink:
    """File wrapper that hashes the bytes actually written, with no ``name`` attribute.

    GzipFile would otherwise copy fileobj.name into the gzip header and make
    the SHA-256 depend on the output path.
    """

    def __init__(self, raw: BinaryIO) -> None:
        self._raw = raw
        self._sha = hashlib.sha256()
        self._closed = False

    def write(self, data: bytes | bytearray | memoryview) -> int:
        self._sha.update(data)
        written = self._raw.write(data)
        return len(data) if written is None else written

    def flush(self) -> None:
        self._raw.flush()

    def close(self) -> None:
        if not self._closed:
            self._raw.close()
            self._closed = True

    def hexdigest(self) -> str:
        return self._sha.hexdigest()


class _Output:
    """Binary FASTQ writer. Gzip members use a fixed mtime and no filename."""

    def __init__(self, path: Path) -> None:
        self._closed = False
        self._gzip = path.name.endswith(".gz")
        self._raw = open(path, "wb")
        self._hasher = _HashingSink(self._raw)
        try:
            if self._gzip:
                self._fh: gzip.GzipFile | _HashingSink = gzip.GzipFile(
                    fileobj=self._hasher,
                    mode="wb",
                    compresslevel=GZIP_COMPRESSLEVEL,
                    mtime=GZIP_MTIME,
                )
            else:
                self._fh = self._hasher
        except Exception:
            self._closed = True
            self._hasher.close()
            raise

    def write_record(self, record: tuple[bytes, bytes, bytes, bytes]) -> None:
        for line in record:
            self._fh.write(line)

    def finish(self) -> str:
        self._close()
        return self._hasher.hexdigest()

    def abort(self) -> None:
        try:
            self._close()
        except Exception:
            self._closed = True

    def _close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self._fh.close()
        finally:
            if self._fh is not self._hasher:
                self._hasher.close()


def emit_pairs(
    r1_path: Path,
    r2_path: Path,
    out1_path: Path,
    out2_path: Path,
    n_pairs: int,
    indices: list[int] | None,
) -> tuple[str, str]:
    """Pass 2: copy selected records. ``indices is None`` keeps every pair.

    Partial outputs are removed if this pass fails.
    """
    cursor = 0
    seen = 0
    out1: _Output | None = None
    out2: _Output | None = None
    completed = False
    try:
        out1 = _Output(out1_path)
        out2 = _Output(out2_path)
        for pair_number, left, right in iter_pair_records(r1_path, r2_path):
            index = pair_number - 1
            keep = indices is None or (cursor < len(indices) and index == indices[cursor])
            if keep:
                out1.write_record(left)
                out2.write_record(right)
                if indices is not None:
                    cursor += 1
            seen += 1
        if seen != n_pairs:
            raise PairSyncError(
                f"pair count changed between passes ({n_pairs} then {seen})"
            )
        if indices is not None and cursor != len(indices):
            raise SubsampleError("internal: did not emit every selected pair")
        digest1 = out1.finish()
        digest2 = out2.finish()
        completed = True
        return digest1, digest2
    finally:
        if not completed:
            if out1 is not None:
                out1.abort()
            if out2 is not None:
                out2.abort()
            out1_path.unlink(missing_ok=True)
            out2_path.unlink(missing_ok=True)


def _require_input(path: Path) -> None:
    if not path.is_file():
        raise SubsampleError(f"input is not a file: {path}")


def _require_output_dir(path: Path) -> None:
    parent = path.parent
    if parent == Path(""):
        parent = Path(".")
    if not parent.is_dir():
        raise SubsampleError(f"output directory does not exist: {parent}")


def _same_path(left: Path, right: Path) -> bool:
    return left.resolve() == right.resolve()


def subsample(args: argparse.Namespace) -> int:
    _require_input(args.r1)
    _require_input(args.r2)
    _require_output_dir(args.out1)
    _require_output_dir(args.out2)

    n_pairs, detected_len, variable = scan_pairs(args.r1, args.r2)
    read_len = detected_len
    clamped = False

    if args.num_pairs is not None:
        mode = "num-pairs"
        requested = args.num_pairs
    elif args.fraction is not None:
        mode = "fraction"
        requested = pairs_for_fraction(args.fraction, n_pairs)
        if requested < 1:
            raise SubsampleError(
                f"--fraction {args.fraction} of {n_pairs} pairs rounds to "
                f"{requested}; at least 1 pair is required"
            )
    else:
        mode = "target-depth"
        if args.read_len is not None:
            read_len = args.read_len
            if args.read_len != detected_len:
                print(
                    f"warning: --read-len {args.read_len} differs from the first "
                    f"R1 read length {detected_len}; the depth formula uses {args.read_len}",
                    file=sys.stderr,
                    flush=True,
                )
        elif variable:
            print(
                f"warning: R1 read lengths are not uniform; the depth formula uses "
                f"the first read ({detected_len} bp)",
                file=sys.stderr,
                flush=True,
            )
        requested = pairs_for_depth(args.target_depth, args.genome_size, read_len)
        if requested < 1:
            raise SubsampleError(
                f"target depth {args.target_depth}x, genome {args.genome_size} bp, "
                f"and read length {read_len} round to {requested} pairs"
            )

    if requested > n_pairs:
        clamped = True
        print(
            f"warning: requested {requested} pairs but only {n_pairs} are available; "
            "keeping all pairs",
            file=sys.stderr,
            flush=True,
        )
        requested = n_pairs

    indices = select_pair_indices(n_pairs, requested, args.seed)
    digest1, digest2 = emit_pairs(
        args.r1,
        args.r2,
        args.out1,
        args.out2,
        n_pairs,
        indices,
    )
    print(
        "[subsample] "
        f"mode={mode} pairs_in={n_pairs} pairs_out={requested} seed={args.seed} "
        f"read_len={read_len} clamped={'yes' if clamped else 'no'} "
        f"r1_sha256={digest1} r2_sha256={digest2}",
        file=sys.stderr,
        flush=True,
    )
    return 0


def _positive_int(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if value < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return value


def _positive_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if value <= 0.0 or value != value or value == float("inf"):
        raise argparse.ArgumentTypeError("must be finite and > 0")
    return value


def _fraction(text: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if value != value or value == float("inf") or not 0.0 < value <= 1.0:
        raise argparse.ArgumentTypeError("must satisfy 0 < fraction <= 1")
    return value


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Deterministically subsample paired FASTQ reads while keeping each "
            "R1/R2 mate pair together. Used for fastBact F26 Stage E coverage gradients."
        ),
        epilog=(
            "examples:\n"
            "  python3 scripts/subsample_fastq_pairs.py \\\n"
            "    --r1 in_R1.fastq.gz --r2 in_R2.fastq.gz \\\n"
            "    --out1 out_R1.fastq.gz --out2 out_R2.fastq.gz \\\n"
            "    --target-depth 50 --genome-size 2800000 --seed 42\n"
            "\n"
            "  python3 scripts/subsample_fastq_pairs.py --self-test\n"
            "\n"
            "k = round((target_depth * genome_size) / (2 * read_len)).\n"
            "If k is greater than the number of input pairs, all pairs are kept."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--r1", type=Path, help="Input R1 FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--r2", type=Path, help="Input R2 FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--out1", type=Path, help="Output R1 path (gzip if it ends with .gz)")
    parser.add_argument("--out2", type=Path, help="Output R2 path (gzip if it ends with .gz)")
    parser.add_argument("--target-depth", type=_positive_float, help="Target coverage depth")
    parser.add_argument("--genome-size", type=_positive_int, help="Genome size in bp")
    parser.add_argument(
        "--read-len",
        type=_positive_int,
        help="Read length for the depth formula (default: first R1 sequence length)",
    )
    parser.add_argument("--num-pairs", type=_positive_int, help="Exact number of pairs to keep")
    parser.add_argument(
        "--fraction",
        type=_fraction,
        help="Fraction of pairs to keep, with 0 < f <= 1",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Deterministic sampling seed (default: {DEFAULT_SEED})",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run the built-in pairing, determinism, and edge-case tests",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.self_test:
        return args

    missing = [name for name in ("r1", "r2", "out1", "out2") if getattr(args, name) is None]
    if missing:
        parser.error(
            "required unless --self-test: " + ", ".join("--" + name for name in missing)
        )

    depth_parts = args.target_depth is not None or args.genome_size is not None
    read_len_set = args.read_len is not None
    modes: list[str] = []
    if depth_parts or read_len_set:
        modes.append("target-depth")
    if args.num_pairs is not None:
        modes.append("num-pairs")
    if args.fraction is not None:
        modes.append("fraction")
    if len(modes) != 1:
        parser.error(
            "specify exactly one of --num-pairs, --fraction, "
            "or (--target-depth with --genome-size)"
        )
    if modes[0] == "target-depth" and (
        args.target_depth is None or args.genome_size is None
    ):
        parser.error("--target-depth and --genome-size must be given together")
    if _same_path(args.out1, args.out2):
        parser.error("--out1 and --out2 must be different paths")
    for source, dest in (
        (args.r1, args.out1),
        (args.r1, args.out2),
        (args.r2, args.out1),
        (args.r2, args.out2),
    ):
        if _same_path(source, dest):
            parser.error(f"refusing to overwrite input {source} with output {dest}")
    return args


def _load_records(path: Path) -> list[tuple[bytes, bytes, bytes, bytes]]:
    records: list[tuple[bytes, bytes, bytes, bytes]] = []
    with open_fastq_input(path) as handle:
        pair_number = 1
        while True:
            record = read_record(handle, path, pair_number)
            if record is None:
                return records
            records.append(record)
            pair_number += 1


def _oracle_indices(n_pairs: int, k: int, seed: int) -> list[int]:
    """Independent copy of the sampling rule used as the test oracle."""
    if k >= n_pairs:
        return list(range(n_pairs))
    return sorted(random.Random(seed).sample(range(n_pairs), k))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _summary_field(text: str, key: str) -> str:
    prefix = key + "="
    for line in text.splitlines():
        if line.startswith("[subsample] "):
            for token in line.split():
                if token.startswith(prefix):
                    return token[len(prefix) :]
    raise AssertionError(f"{key} missing from subsample summary: {text!r}")


def _seq(index: int, mate: int, read_len: int) -> bytes:
    token = f"{index:08d}".encode("ascii")
    if read_len < len(token):
        raise RuntimeError("fixture read length is shorter than the index token")
    pad = b"A" if mate == 1 else b"C"
    return pad * (read_len - len(token)) + token


def _headers(index: int, style: str) -> tuple[bytes, bytes]:
    if style == "slash":
        return f"@pair{index:06d}/1".encode(), f"@pair{index:06d}/2".encode()
    if style == "illumina":
        base = f"@INST:1:FC:1:1:{index}:1234"
        return f"{base} 1:N:0:ACGT".encode(), f"{base} 2:N:0:ACGT".encode()
    if style == "identical":
        header = f"@spot{index:06d}".encode()
        return header, header
    if style == "both":
        return (
            f"@pair{index:06d}/1 1:N:0:ACGT".encode(),
            f"@pair{index:06d}/2 2:N:0:ACGT".encode(),
        )
    if style == "oneside":
        return f"@pair{index:06d}/1".encode(), f"@pair{index:06d}".encode()
    if style == "conflict":
        return (
            f"@pair{index:06d}/1 2:N:0:ACGT".encode(),
            f"@pair{index:06d}/2 2:N:0:ACGT".encode(),
        )
    if style == "swap":
        base = f"@INST:1:FC:1:1:{index}:1234"
        return f"{base} 2:N:0:ACGT".encode(), f"{base} 1:N:0:ACGT".encode()
    if style == "mismatch":
        return f"@alpha{index:06d}/1".encode(), f"@beta{index:06d}/2".encode()
    raise AssertionError(f"unknown header style {style}")


def _record(header: bytes, seq: bytes, newline: bytes = b"\n") -> bytes:
    return b"".join((header, newline, seq, newline, b"+", newline, b"I" * len(seq), newline))


def _build_pairs(
    n_pairs: int,
    style: str = "slash",
    read_len: int = 20,
    lengths: Sequence[int] | None = None,
) -> tuple[list[bytes], list[bytes]]:
    r1_records: list[bytes] = []
    r2_records: list[bytes] = []
    for index in range(n_pairs):
        this_len = read_len if lengths is None else lengths[index]
        header1, header2 = _headers(index, style)
        r1_records.append(_record(header1, _seq(index, 1, this_len)))
        r2_records.append(_record(header2, _seq(index, 2, this_len)))
    return r1_records, r2_records


def _write_bytes(path: Path, payload: bytes) -> None:
    if path.name.endswith(".gz"):
        with gzip.GzipFile(filename=str(path), mode="wb", compresslevel=6, mtime=0) as handle:
            handle.write(payload)
    else:
        path.write_bytes(payload)


def _write_pairs(path1: Path, path2: Path, r1_records: Sequence[bytes], r2_records: Sequence[bytes]) -> None:
    _write_bytes(path1, b"".join(r1_records))
    _write_bytes(path2, b"".join(r2_records))


def self_test() -> int:
    failures: list[str] = []

    def check(condition: bool, message: str) -> None:
        if condition:
            print(f"ok   {message}")
        else:
            failures.append(message)
            print(f"FAIL {message}")

    def invoke(argv: Sequence[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        try:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                code = main(list(argv))
        except SystemExit as exc:
            code = 0 if exc.code is None else int(exc.code)
        return code, stdout.getvalue(), stderr.getvalue()

    check(pairs_for_depth(30.5, 5_000_000, 150) == 508_333, "depth formula 30.5x / 5 Mb / 150 bp")
    check(pairs_for_depth(1, 3, 1) == 2, "round half to even: 1.5 -> 2")
    check(pairs_for_depth(1, 1, 1) == 0, "round half to even: 0.5 -> 0")
    check(pairs_for_fraction(0.25, 32) == 8, "fraction 0.25 of 32 pairs is 8")
    oracle_8 = _oracle_indices(32, 8, 42)
    check(
        select_pair_indices(32, 8, 42) == oracle_8,
        "select_pair_indices matches sorted Random(42).sample(range(32), 8)",
    )
    check(select_pair_indices(32, 32, 42) is None, "selecting every pair does not build a permutation")
    check(oracle_8 != list(range(8)), "seed 42 is not a prefix sample")
    check(
        select_pair_indices(32, 8, -1) != select_pair_indices(32, 8, 42),
        "a negative seed is accepted and changes the sample",
    )
    for seed in (0, 1, 42, -3, 99):
        identity = sorted(random.Random(seed).sample(range(16), 16))
        check(identity == list(range(16)), f"sorted full sample is the original order for seed {seed}")

    with tempfile.TemporaryDirectory(prefix="f26-subsample-self-test-") as tmp_name:
        tmp = Path(tmp_name)
        n_pairs = 32
        r1_records, r2_records = _build_pairs(n_pairs, "slash", 20)
        r1_path = tmp / "in_R1.fastq"
        r2_path = tmp / "in_R2.fastq"
        _write_pairs(r1_path, r2_path, r1_records, r2_records)

        def run_case(
            name: str,
            extra: Sequence[str],
            out1: Path,
            out2: Path,
            in1: Path = r1_path,
            in2: Path = r2_path,
            seed: int | None = 42,
        ) -> tuple[int, str, str]:
            argv = ["--r1", str(in1), "--r2", str(in2), "--out1", str(out1), "--out2", str(out2)]
            if seed is not None:
                argv.extend(["--seed", str(seed)])
            argv.extend(extra)
            code, stdout, stderr = invoke(argv)
            if code != 0:
                check(False, f"{name} exit {code}: {stderr.strip()}")
            return code, stdout, stderr

        out1 = tmp / "a_R1.fastq"
        out2 = tmp / "a_R2.fastq"
        code, stdout, stderr = run_case("num-pairs", ["--num-pairs", "8"], out1, out2)
        check(code == 0 and stdout == "", "num-pairs exits 0 and writes nothing to stdout")
        kept1 = _load_records(out1)
        kept2 = _load_records(out2)
        expected = [r1_records[i] for i in oracle_8]
        expected2 = [r2_records[i] for i in oracle_8]
        check(kept1 == [tuple(_record_parts(item)) for item in expected], "R1 records are the sampled input pairs in file order")
        check(kept2 == [tuple(_record_parts(item)) for item in expected2], "R2 records stay with the same sampled indices")
        check(
            all(left[0].split()[0][:-2] == right[0].split()[0][:-2] for left, right in zip(kept1, kept2)),
            "every emitted pair still shares its read id after /1 and /2 removal",
        )
        check(_summary_field(stderr, "r1_sha256") == _sha256_file(out1), "reported R1 SHA-256 matches the output file")
        check(_summary_field(stderr, "r2_sha256") == _sha256_file(out2), "reported R2 SHA-256 matches the output file")
        check(_sha256_file(out1) != _sha256_file(out2), "R1 and R2 outputs are not copies of each other")

        code_again, _stdout_again, stderr_again = run_case("rerun", ["--num-pairs", "8"], out1, out2)
        check(
            code_again == 0 and _sha256_file(out1) == _summary_field(stderr, "r1_sha256"),
            "re-running with the same seed keeps the R1 SHA-256",
        )
        check(
            _summary_field(stderr_again, "r1_sha256") == _summary_field(stderr, "r1_sha256")
            and _summary_field(stderr_again, "r2_sha256") == _summary_field(stderr, "r2_sha256"),
            "re-running with the same seed produces identical R1 and R2 SHA-256",
        )

        other1 = tmp / "other" / "b_R1.fastq.gz"
        other2 = tmp / "other" / "b_R2.fastq.gz"
        other1.parent.mkdir()
        gz_code, _gz_stdout, gz_stderr = run_case(
            "gzip-path",
            ["--num-pairs", "8"],
            other1,
            other2,
        )
        plain_sha_1 = _sha256_file(out1)
        check(gz_code == 0, "gzip outputs are written when the path ends in .gz")
        check(
            _load_records(other1) == kept1 and _load_records(other2) == kept2,
            "gzip outputs decompress to the same pairs as plain outputs",
        )
        third1 = tmp / "elsewhere_R1.fastq.gz"
        third2 = tmp / "elsewhere_R2.fastq.gz"
        third_code, _third_stdout, third_stderr = run_case(
            "gzip-elsewhere",
            ["--num-pairs", "8"],
            third1,
            third2,
        )
        check(
            third_code == 0
            and _sha256_file(third1) == _sha256_file(other1)
            and _sha256_file(third2) == _sha256_file(other2)
            and _summary_field(gz_stderr, "r1_sha256") == _summary_field(third_stderr, "r1_sha256"),
            "gzip SHA-256 is stable across output paths",
        )
        check(plain_sha_1 != _sha256_file(other1), "plain and gzip containers differ, as expected")

        div1 = tmp / "div_R1.fastq"
        div2 = tmp / "div_R2.fastq"
        div_code, _div_stdout, div_stderr = run_case(
            "seed-43",
            ["--num-pairs", "8"],
            div1,
            div2,
            seed=43,
        )
        check(
            div_code == 0
            and _summary_field(div_stderr, "r1_sha256") != _summary_field(stderr, "r1_sha256")
            and _summary_field(div_stderr, "r2_sha256") != _summary_field(stderr, "r2_sha256"),
            "a different seed produces different R1 and R2 SHA-256",
        )
        check(
            _load_records(div1) != kept1,
            "a different seed keeps a different set of pairs",
        )

        def_code, _def_stdout, def_stderr = run_case(
            "default-seed",
            ["--num-pairs", "8"],
            tmp / "def_R1.fastq",
            tmp / "def_R2.fastq",
            seed=None,
        )
        check(
            def_code == 0
            and _summary_field(def_stderr, "seed") == str(DEFAULT_SEED)
            and _summary_field(def_stderr, "r1_sha256") == _summary_field(stderr, "r1_sha256"),
            "omitted --seed uses 42 and matches an explicit seed of 42",
        )

        all1 = tmp / "all_R1.fastq"
        all2 = tmp / "all_R2.fastq"
        all_code, _all_stdout, all_stderr = run_case(
            "num-above",
            ["--num-pairs", "10000"],
            all1,
            all2,
            seed=7,
        )
        check(all_code == 0, "target pair count above the input exits 0")
        check("keeping all pairs" in all_stderr, "target above the input reports that all pairs are kept")
        check(
            _load_records(all1) == [tuple(_record_parts(item)) for item in r1_records]
            and _load_records(all2) == [tuple(_record_parts(item)) for item in r2_records],
            "target above the input keeps every pair in the original order",
        )
        all_b1 = tmp / "allb_R1.fastq"
        all_b2 = tmp / "allb_R2.fastq"
        all_b_code, _all_b_stdout, all_b_stderr = run_case(
            "num-above-other-seed",
            ["--num-pairs", "10000"],
            all_b1,
            all_b2,
            seed=99,
        )
        check(
            all_b_code == 0
            and _sha256_file(all_b1) == _sha256_file(all1)
            and _summary_field(all_b_stderr, "pairs_out") == str(n_pairs),
            "keeping all pairs does not depend on the seed",
        )

        eq1 = tmp / "eq_R1.fastq"
        eq2 = tmp / "eq_R2.fastq"
        eq_code, _eq_stdout, eq_stderr = run_case(
            "num-equal",
            ["--num-pairs", str(n_pairs)],
            eq1,
            eq2,
        )
        check(
            eq_code == 0
            and _summary_field(eq_stderr, "clamped") == "no"
            and _load_records(eq1) == [tuple(_record_parts(item)) for item in r1_records],
            "requesting exactly N pairs keeps the file and is not reported as clamped",
        )

        depth1 = tmp / "depth_R1.fastq"
        depth2 = tmp / "depth_R2.fastq"
        # k = round((2 * 200) / (2 * 20)) = 10
        depth_code, _depth_stdout, depth_stderr = run_case(
            "depth-auto",
            ["--target-depth", "2", "--genome-size", "200"],
            depth1,
            depth2,
        )
        depth_oracle = _oracle_indices(n_pairs, 10, 42)
        check(
            depth_code == 0
            and _summary_field(depth_stderr, "read_len") == "20"
            and _summary_field(depth_stderr, "pairs_out") == "10"
            and _load_records(depth1) == [tuple(_record_parts(r1_records[i])) for i in depth_oracle],
            "target depth detects read length from the first R1 read",
        )

        forced1 = tmp / "forced_R1.fastq"
        forced2 = tmp / "forced_R2.fastq"
        # explicit length 10: k = round(400 / 20) = 20
        forced_code, _forced_stdout, forced_stderr = run_case(
            "depth-explicit",
            ["--target-depth", "2", "--genome-size", "200", "--read-len", "10"],
            forced1,
            forced2,
        )
        check(
            forced_code == 0
            and "differs from the first R1 read length" in forced_stderr
            and _summary_field(forced_stderr, "pairs_out") == "20"
            and _summary_field(forced_stderr, "read_len") == "10",
            "--read-len overrides detection and is what the depth formula uses",
        )

        sat1 = tmp / "sat_R1.fastq.gz"
        sat2 = tmp / "sat_R2.fastq.gz"
        sat_code, _sat_stdout, sat_stderr = run_case(
            "depth-saturated",
            ["--target-depth", "150", "--genome-size", "1000000", "--read-len", "20"],
            sat1,
            sat2,
        )
        check(
            sat_code == 0
            and _summary_field(sat_stderr, "clamped") == "yes"
            and _summary_field(sat_stderr, "pairs_out") == str(n_pairs)
            and _load_records(sat1) == [tuple(_record_parts(item)) for item in r1_records],
            "a 150x target above the available reads keeps every pair and exits 0",
        )

        frac1 = tmp / "frac_R1.fastq"
        frac2 = tmp / "frac_R2.fastq"
        frac_code, _frac_stdout, frac_stderr = run_case(
            "fraction",
            ["--fraction", "0.25"],
            frac1,
            frac2,
        )
        check(
            frac_code == 0
            and _load_records(frac1) == kept1
            and _summary_field(frac_stderr, "mode") == "fraction",
            "fraction 0.25 matches the 8-pair sample for this file",
        )
        full1 = tmp / "full_R1.fastq"
        full2 = tmp / "full_R2.fastq"
        full_code, _full_stdout, full_stderr = run_case(
            "fraction-one",
            ["--fraction", "1"],
            full1,
            full2,
        )
        check(
            full_code == 0 and _summary_field(full_stderr, "pairs_out") == str(n_pairs),
            "fraction 1 keeps every pair",
        )

        var_r1, var_r2 = _build_pairs(8, "slash", lengths=[20, 10, 10, 10, 10, 10, 10, 10])
        var_in1 = tmp / "var_R1.fastq"
        var_in2 = tmp / "var_R2.fastq"
        _write_pairs(var_in1, var_in2, var_r1, var_r2)
        var_code, _var_stdout, var_stderr = run_case(
            "variable",
            ["--target-depth", "2", "--genome-size", "100"],
            tmp / "var_out_R1.fastq",
            tmp / "var_out_R2.fastq",
            in1=var_in1,
            in2=var_in2,
        )
        # first read is 20 bp: k = round(200 / 40) = 5, not the shorter length
        check(
            var_code == 0
            and "not uniform" in var_stderr
            and _summary_field(var_stderr, "pairs_out") == "5"
            and _summary_field(var_stderr, "read_len") == "20",
            "variable-length input uses the first R1 read for the depth formula",
        )

        for style in ("illumina", "identical", "both", "oneside"):
            style_r1, style_r2 = _build_pairs(6, style, 20)
            style_in1 = tmp / f"{style}_R1.fastq"
            style_in2 = tmp / f"{style}_R2.fastq"
            _write_pairs(style_in1, style_in2, style_r1, style_r2)
            style_out1 = tmp / f"{style}_out_R1.fastq"
            style_out2 = tmp / f"{style}_out_R2.fastq"
            style_code, _style_stdout, _style_stderr = run_case(
                style,
                ["--num-pairs", "3"],
                style_out1,
                style_out2,
                in1=style_in1,
                in2=style_in2,
            )
            style_oracle = _oracle_indices(6, 3, 42)
            parsed1 = _load_records(style_out1)
            parsed2 = _load_records(style_out2)
            ids_ok = True
            for record1, record2 in zip(parsed1, parsed2):
                try:
                    assert_pair_headers(record1[0], record2[0], 1)
                except SubsampleError:
                    ids_ok = False
            check(
                style_code == 0
                and ids_ok
                and parsed1 == [tuple(_record_parts(style_r1[i])) for i in style_oracle]
                and parsed2 == [tuple(_record_parts(style_r2[i])) for i in style_oracle],
                f"{style} headers stay paired and the sampled records are unchanged",
            )

        gz_in1 = tmp / "gz_in_R1.fastq.gz"
        gz_in2 = tmp / "gz_in_R2.fastq.gz"
        _write_pairs(gz_in1, gz_in2, r1_records, r2_records)
        mixed_code, _mixed_stdout, _mixed_stderr = run_case(
            "gzip-input",
            ["--num-pairs", "8"],
            tmp / "from_gz_R1.fastq",
            tmp / "from_gz_R2.fastq.gz",
            in1=gz_in1,
            in2=gz_in2,
        )
        check(
            mixed_code == 0
            and _load_records(tmp / "from_gz_R1.fastq") == kept1
            and _load_records(tmp / "from_gz_R2.fastq.gz") == kept2,
            "gzip inputs and mixed output containers preserve the sampled pairs",
        )

        crlf1 = tmp / "crlf_R1.fastq"
        crlf2 = tmp / "crlf_R2.fastq"
        crlf_r1 = _record(b"@read/1", b"ACGTACGTACGTACGTACGT", b"\r\n")
        crlf_r2 = _record(b"@read/2", b"TGCATGCATGCATGCATGCA", b"\r\n")
        crlf1.write_bytes(crlf_r1)
        crlf2.write_bytes(crlf_r2)
        crlf_out1 = tmp / "crlf_out_R1.fastq"
        crlf_out2 = tmp / "crlf_out_R2.fastq"
        crlf_code, _crlf_stdout, _crlf_stderr = run_case(
            "crlf",
            ["--num-pairs", "1"],
            crlf_out1,
            crlf_out2,
            in1=crlf1,
            in2=crlf2,
        )
        check(
            crlf_code == 0
            and crlf_out1.read_bytes() == crlf_r1
            and crlf_out2.read_bytes() == crlf_r2,
            "CRLF records are copied unchanged and /1 /2 headers still match",
        )

        def expect_fail(name: str, argv: Sequence[str], expected: int, sentinel: Path | None = None) -> None:
            before = sentinel.read_bytes() if sentinel is not None else None
            code, _stdout, stderr = invoke(argv)
            intact = True if sentinel is None else sentinel.read_bytes() == before
            check(code == expected and intact, f"{name} exits {expected} (got {code}); {stderr.strip()}")

        sentinel = tmp / "sentinel_R1.fastq"
        sentinel.write_bytes(b"KEEP\n")
        other_sentinel = tmp / "sentinel_R2.fastq"
        other_sentinel.write_bytes(b"KEEP\n")
        bad_r1, bad_r2 = _build_pairs(2, "mismatch")
        bad_in1 = tmp / "bad_R1.fastq"
        bad_in2 = tmp / "bad_R2.fastq"
        _write_pairs(bad_in1, bad_in2, bad_r1, bad_r2)
        expect_fail(
            "mismatched read ids",
            ["--r1", str(bad_in1), "--r2", str(bad_in2), "--out1", str(sentinel), "--out2", str(other_sentinel), "--num-pairs", "1"],
            1,
            sentinel,
        )
        check(other_sentinel.read_bytes() == b"KEEP\n", "a header mismatch does not replace an existing R2 output")

        for style, label in (("swap", "swapped Illumina mate tags"), ("conflict", "conflicting /1 and Illumina tag")):
            bad1, bad2 = _build_pairs(2, style)
            in1 = tmp / f"{style}_in_R1.fastq"
            in2 = tmp / f"{style}_in_R2.fastq"
            _write_pairs(in1, in2, bad1, bad2)
            expect_fail(
                label,
                ["--r1", str(in1), "--r2", str(in2), "--out1", str(tmp / f"{style}_o1.fastq"), "--out2", str(tmp / f"{style}_o2.fastq"), "--num-pairs", "1"],
                1,
            )
            check(not (tmp / f"{style}_o1.fastq").exists(), f"{label} does not leave an output file")

        short_r1, short_r2 = _build_pairs(3, "slash")
        _write_pairs(tmp / "uneq_R1.fastq", tmp / "uneq_R2.fastq", short_r1, short_r2[:-1])
        expect_fail(
            "unequal pair counts",
            ["--r1", str(tmp / "uneq_R1.fastq"), "--r2", str(tmp / "uneq_R2.fastq"), "--out1", str(tmp / "uneq_o1.fastq"), "--out2", str(tmp / "uneq_o2.fastq"), "--num-pairs", "1"],
            1,
        )

        (tmp / "empty_R1.fastq").write_bytes(b"")
        (tmp / "empty_R2.fastq").write_bytes(b"")
        expect_fail(
            "empty inputs",
            ["--r1", str(tmp / "empty_R1.fastq"), "--r2", str(tmp / "empty_R2.fastq"), "--out1", str(tmp / "empty_o1.fastq"), "--out2", str(tmp / "empty_o2.fastq"), "--num-pairs", "1"],
            1,
        )
        (tmp / "trunc_R1.fastq").write_bytes(b"@read/1\nACGT\n+\n")
        (tmp / "trunc_R2.fastq").write_bytes(_record(b"@read/2", b"ACGT"))
        expect_fail(
            "truncated R1",
            ["--r1", str(tmp / "trunc_R1.fastq"), "--r2", str(tmp / "trunc_R2.fastq"), "--out1", str(tmp / "trunc_o1.fastq"), "--out2", str(tmp / "trunc_o2.fastq"), "--num-pairs", "1"],
            1,
        )
        (tmp / "badqual_R1.fastq").write_bytes(b"@read/1\nACGT\n+\nII\n")
        (tmp / "badqual_R2.fastq").write_bytes(_record(b"@read/2", b"ACGT"))
        expect_fail(
            "sequence/quality length mismatch",
            ["--r1", str(tmp / "badqual_R1.fastq"), "--r2", str(tmp / "badqual_R2.fastq"), "--out1", str(tmp / "bq_o1.fastq"), "--out2", str(tmp / "bq_o2.fastq"), "--num-pairs", "1"],
            1,
        )
        (tmp / "badplus_R1.fastq").write_bytes(b"@read/1\nACGT\n-\nIIII\n")
        (tmp / "badplus_R2.fastq").write_bytes(_record(b"@read/2", b"ACGT"))
        expect_fail(
            "missing plus separator",
            ["--r1", str(tmp / "badplus_R1.fastq"), "--r2", str(tmp / "badplus_R2.fastq"), "--out1", str(tmp / "bp_o1.fastq"), "--out2", str(tmp / "bp_o2.fastq"), "--num-pairs", "1"],
            1,
        )
        (tmp / "not_R1.fastq.gz").write_bytes(b"this is not gzip")
        (tmp / "not_R2.fastq.gz").write_bytes(b"this is not gzip")
        expect_fail(
            "invalid gzip",
            ["--r1", str(tmp / "not_R1.fastq.gz"), "--r2", str(tmp / "not_R2.fastq.gz"), "--out1", str(tmp / "ng_o1.fastq"), "--out2", str(tmp / "ng_o2.fastq"), "--num-pairs", "1"],
            1,
        )
        expect_fail(
            "missing output directory",
            ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "missing-dir" / "o1.fastq"), "--out2", str(tmp / "o2.fastq"), "--num-pairs", "1"],
            1,
        )
        tiny_r1, tiny_r2 = _build_pairs(3, "slash", 20)
        _write_pairs(tmp / "tiny_R1.fastq", tmp / "tiny_R2.fastq", tiny_r1, tiny_r2)
        zero_code, _zero_stdout, zero_stderr = invoke(
            [
                "--r1", str(tmp / "tiny_R1.fastq"),
                "--r2", str(tmp / "tiny_R2.fastq"),
                "--out1", str(tmp / "zero_o1.fastq"),
                "--out2", str(tmp / "zero_o2.fastq"),
                "--fraction", "0.1",
            ]
        )
        check(
            zero_code == 1 and "rounds to 0" in zero_stderr and not (tmp / "zero_o1.fastq").exists(),
            "a fraction that rounds to 0 pairs is rejected",
        )
        depth0_code, _depth0_stdout, depth0_stderr = invoke(
            [
                "--r1", str(tmp / "tiny_R1.fastq"),
                "--r2", str(tmp / "tiny_R2.fastq"),
                "--out1", str(tmp / "d0_o1.fastq"),
                "--out2", str(tmp / "d0_o2.fastq"),
                "--target-depth", "1",
                "--genome-size", "1",
                "--read-len", "1",
            ]
        )
        check(
            depth0_code == 1 and "round to 0" in depth0_stderr,
            "a depth that rounds to 0 pairs is rejected",
        )

        usage_cases = [
            ("missing required arguments", ["--num-pairs", "1"]),
            ("two subsampling modes", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--num-pairs", "1", "--fraction", "0.5"]),
            ("fraction 0", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--fraction", "0"]),
            ("fraction above 1", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--fraction", "1.5"]),
            ("depth without genome size", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--target-depth", "30"]),
            ("read-len without depth", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--read-len", "150", "--num-pairs", "1"]),
            ("identical outputs", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "same.fastq"), "--out2", str(tmp / "same.fastq"), "--num-pairs", "1"]),
            ("output overwrites input", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(r1_path), "--out2", str(tmp / "u2.fastq"), "--num-pairs", "1"]),
            ("num-pairs of 0", ["--r1", str(r1_path), "--r2", str(r2_path), "--out1", str(tmp / "u1.fastq"), "--out2", str(tmp / "u2.fastq"), "--num-pairs", "0"]),
        ]
        for label, argv in usage_cases:
            code, _stdout, _stderr = invoke(argv)
            check(code == 2, f"{label} exits 2 (got {code})")

    if failures:
        print(f"SELF_TEST_FAIL {len(failures)}")
        return 1
    print("SELF_TEST_PASS")
    return 0


def _record_parts(payload: bytes) -> tuple[bytes, bytes, bytes, bytes]:
    lines = payload.splitlines(keepends=True)
    if len(lines) != 4:
        raise AssertionError(f"fixture record does not have 4 lines: {payload!r}")
    return lines[0], lines[1], lines[2], lines[3]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        return self_test()
    try:
        return subsample(args)
    except SubsampleError as exc:
        print(f"error: {exc}", file=sys.stderr, flush=True)
        return 1
    except (OSError, EOFError, zlib.error) as exc:
        print(f"error: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
