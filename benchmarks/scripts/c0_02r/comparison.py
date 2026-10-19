from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import os
import shutil
import tempfile
from typing import TypedDict
from .alignment import (
    Alignment,
    AlignmentError,
    Event,
    ambiguous,
    parse_paf,
    project_query,
)
from .execution import Executor, require_compute
from .structure import extract_events


def minimap2() -> str:
    candidate = os.environ.get("MINIMAP2") or shutil.which("minimap2")
    if candidate:
        return candidate
    bundled = (
        Path.home() / ".local/lib/python3.9/site-packages/quast_libs/minimap2/minimap2"
    )
    if bundled.is_file():
        return str(bundled)
    raise AlignmentError("Set MINIMAP2 to an executable minimap2 path")


def align(
    executor: Executor, pair: tuple[Path, Path], label: str
) -> tuple[Alignment, ...]:
    target, query = pair
    result = executor.run(
        (
            minimap2(),
            "-x",
            "asm5",
            "-c",
            "--cs=long",
            "--secondary=yes",
            "-N",
            "20",
            "-t",
            "8",
            str(target),
            str(query),
        ),
        label,
    )
    return parse_paf(Path(result.stdout).read_text())


def fasta_lengths(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    name = ""
    with path.open() as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in sizes:
                    raise AlignmentError(f"Duplicate FASTA name: {name}")
                sizes[name] = 0
            elif line:
                if not name:
                    raise AlignmentError(f"FASTA sequence before header in {path}")
                sizes[name] += len(line)
    if not sizes or any(length == 0 for length in sizes.values()):
        raise AlignmentError(f"Empty FASTA: {path}")
    return sizes


def unaligned_events(
    pair: tuple[Path, Path], alignments: tuple[Alignment, ...]
) -> list[Event]:
    events: list[Event] = []
    for side, path in enumerate(pair):
        for name, length in fasta_lengths(path).items():
            spans = sorted(
                (a.target_start, a.target_end)
                if side == 0
                else (a.query_start, a.query_end)
                for a in alignments
                if a.primary and (a.target if side == 0 else a.query) == name
            )
            end = 0
            for start, stop in [*spans, (length, length)]:
                if start > end:
                    events.append(
                        Event(
                            "UNALIGNED_BASELINE_REGION"
                            if side == 0
                            else "UNALIGNED_FIXED_REGION",
                            name if side == 0 else "NA",
                            end if side == 0 else 0,
                            start if side == 0 else 0,
                            name if side == 1 else "NA",
                            end if side == 1 else 0,
                            start if side == 1 else 0,
                            ".",
                            start - end,
                            evidence="NO_PRIMARY_ALIGNMENT",
                        )
                    )
                end = max(end, stop)
    return events


def annotate_reference(event: Event, refs: tuple[Alignment, ...]) -> Event:
    if event.evidence == "AMBIGUOUS_REPETITIVE_ALIGNMENT":
        return event
    hits: list[tuple[str, int]] = []
    for a in refs:
        if a.primary and a.query == event.fixed_contig and not ambiguous(a, refs):
            if event.fixed_end > event.fixed_start and event.fixed_end > a.query_end:
                continue
            point = project_query(a, event.fixed_start)
            if point is not None:
                if event.fixed_start == event.fixed_end and a.strand == "-":
                    point += 1
                hits.append((a.target, point))
    if len(hits) == 1:
        name, point = hits[0]
        return replace(
            event,
            ref_coords=f"{name}:{point}",
            evidence=event.evidence + ";REFERENCE_COORDINATE_ONLY",
        )
    return event


class LegacyEvent(TypedDict):
    event_type: str
    event_len: int
    baseline_contig: str
    baseline_start: int
    baseline_end: int
    fixed_contig: str
    fixed_start: int
    fixed_end: int
    evidence: str


def compare_fasta(old_p: str, new_p: str) -> list[LegacyEvent]:
    require_compute()
    with tempfile.TemporaryDirectory(prefix="c002r-alignment-") as temporary:
        pair = (Path(old_p), Path(new_p))
        alignments = align(Executor(Path(temporary), []), pair, "pair")
        events = extract_events(alignments) + unaligned_events(pair, alignments)
        return [
            {
                "event_type": event.event_type,
                "event_len": event.event_len,
                "baseline_contig": event.baseline_contig,
                "baseline_start": event.baseline_start,
                "baseline_end": event.baseline_end,
                "fixed_contig": event.fixed_contig,
                "fixed_start": event.fixed_start,
                "fixed_end": event.fixed_end,
                "evidence": event.evidence,
            }
            for event in events
        ]
