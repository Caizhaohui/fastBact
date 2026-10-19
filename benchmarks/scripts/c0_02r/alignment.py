from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import re
from typing import Final, assert_never


class AlignmentError(RuntimeError):
    pass


class Operation(StrEnum):
    MATCH = "match"
    SUBSTITUTION = "SUBSTITUTION"
    INSERTION = "INSERTION"
    DELETION = "DELETION"


@dataclass(frozen=True, slots=True)
class Edit:
    kind: Operation
    target: str
    query: str
    length: int


@dataclass(frozen=True, slots=True)
class Alignment:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int
    mapq: int
    primary: bool
    edits: tuple[Edit, ...]


@dataclass(frozen=True, slots=True)
class Event:
    event_type: str
    baseline_contig: str
    baseline_start: int
    baseline_end: int
    fixed_contig: str
    fixed_start: int
    fixed_end: int
    strand: str
    event_len: int
    baseline_allele: str = ""
    fixed_allele: str = ""
    evidence: str = "PAIRWISE_ALIGNMENT"
    ref_coords: str = "NA"


TOKEN: Final = re.compile(
    r":[0-9]+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+"
)


def parse_edits(cs: str) -> tuple[Edit, ...]:
    tokens = TOKEN.findall(cs)
    if "".join(tokens) != cs or not tokens:
        raise AlignmentError(f"Unsupported or malformed cs tag: {cs[:80]}")
    edits: list[Edit] = []
    for token in tokens:
        match token[0]:
            case ":":
                edits.append(Edit(Operation.MATCH, "", "", int(token[1:])))
            case "=":
                edits.append(
                    Edit(
                        Operation.MATCH,
                        token[1:].upper(),
                        token[1:].upper(),
                        len(token) - 1,
                    )
                )
            case "*":
                edits.append(
                    Edit(Operation.SUBSTITUTION, token[1].upper(), token[2].upper(), 1)
                )
            case "+":
                edits.append(
                    Edit(Operation.INSERTION, "", token[1:].upper(), len(token) - 1)
                )
            case "-":
                edits.append(
                    Edit(Operation.DELETION, token[1:].upper(), "", len(token) - 1)
                )
            case invalid:
                raise AlignmentError(f"Invalid cs operation: {invalid}")
    return tuple(edits)


def parse_paf(text: str) -> tuple[Alignment, ...]:
    alignments: list[Alignment] = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            raise AlignmentError("PAF record has fewer than 12 fields")
        tags = {field[:2]: field[5:] for field in fields[12:]}
        if fields[4] not in ("+", "-") or "cs" not in tags:
            raise AlignmentError("PAF needs strand and cs tag")
        alignment_type = tags.get("tp", "P")
        if alignment_type not in {"P", "I", "S", "i"}:
            raise AlignmentError(f"Unknown minimap2 alignment type: {alignment_type}")
        a = Alignment(
            fields[0],
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            fields[4],
            fields[5],
            int(fields[6]),
            int(fields[7]),
            int(fields[8]),
            int(fields[11]),
            alignment_type in {"P", "I"},
            parse_edits(tags["cs"]),
        )
        if not (
            0 <= a.query_start < a.query_end <= a.query_length
            and 0 <= a.target_start < a.target_end <= a.target_length
            and 0 <= a.mapq <= 255
        ):
            raise AlignmentError("PAF coordinates or mapq outside bounds")
        query_span = sum(
            edit.length for edit in a.edits if edit.kind != Operation.DELETION
        )
        target_span = sum(
            edit.length for edit in a.edits if edit.kind != Operation.INSERTION
        )
        if (
            query_span != a.query_end - a.query_start
            or target_span != a.target_end - a.target_start
        ):
            raise AlignmentError("cs consumption disagrees with PAF coordinates")
        alignments.append(a)
    return tuple(alignments)


def ambiguous(a: Alignment, alignments: tuple[Alignment, ...]) -> bool:
    if a.mapq < 20 or a.mapq == 255:
        return True
    return any(
        b != a
        and b.query == a.query
        and (b.target, b.target_start, b.strand) != (a.target, a.target_start, a.strand)
        and min(a.query_end, b.query_end) - max(a.query_start, b.query_start)
        >= (a.query_end - a.query_start) * 0.5
        for b in alignments
    )


def alignment_events(a: Alignment, evidence: str) -> list[Event]:
    target, query = a.target_start, a.query_start if a.strand == "+" else a.query_end
    direction = 1 if a.strand == "+" else -1
    events: list[Event] = []
    if a.strand == "-":
        events.append(
            Event(
                "REVERSE_COMPLEMENT_ALIGNMENT",
                a.target,
                a.target_start,
                a.target_end,
                a.query,
                a.query_start,
                a.query_end,
                a.strand,
                a.query_end - a.query_start,
                evidence=evidence,
            )
        )
    for edit in a.edits:
        match edit.kind:
            case Operation.MATCH:
                target += edit.length
                query += direction * edit.length
                continue
            case Operation.SUBSTITUTION:
                dt, dq = 1, 1
            case Operation.INSERTION:
                dt, dq = 0, edit.length
            case Operation.DELETION:
                dt, dq = edit.length, 0
            case unreachable:
                assert_never(unreachable)
        next_query = query + direction * dq
        events.append(
            Event(
                edit.kind.value,
                a.target,
                target,
                target + dt,
                a.query,
                min(query, next_query),
                max(query, next_query),
                a.strand,
                edit.length,
                edit.target,
                edit.query,
                evidence,
            )
        )
        target += dt
        query = next_query
    return events


def project_query(a: Alignment, position: int) -> int | None:
    target, query = a.target_start, a.query_start if a.strand == "+" else a.query_end
    direction = 1 if a.strand == "+" else -1
    for edit in a.edits:
        match edit.kind:
            case Operation.MATCH:
                low, high = sorted((query, query + direction * edit.length))
                if low <= position < high:
                    offset = (
                        position - query if direction == 1 else query - 1 - position
                    )
                    return target + offset
                target += edit.length
                query += direction * edit.length
            case Operation.SUBSTITUTION:
                target += 1
                query += direction
            case Operation.INSERTION:
                query += direction * edit.length
            case Operation.DELETION:
                target += edit.length
            case unreachable:
                assert_never(unreachable)
    return None
