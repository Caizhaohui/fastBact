from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Final
from .execution import GateError


QUAST_FIELDS: Final = (
    "N50",
    "NGA50",
    "Genome fraction (%)",
    "# misassemblies",
    "# local misassemblies",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp",
)


@dataclass(frozen=True, slots=True)
class Metrics:
    n50: float
    nga50: float
    genome_fraction: float
    extensive_misassemblies: float
    local_misassemblies: float
    mismatches_per_100kbp: float
    indels_per_100kbp: float


def parse_quast(report: Path) -> Metrics:
    fields: dict[str, str] = {}
    for line in report.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) != 2 or parts[0] in fields:
            raise GateError(
                f"Malformed or duplicate QUAST row in {report}: {line[:100]}"
            )
        fields[parts[0]] = parts[1]
    numbers: list[float] = []
    for name in QUAST_FIELDS:
        if name not in fields:
            raise GateError(f"Missing QUAST metric {name} in {report}")
        try:
            number = float(fields[name])
        except ValueError as error:
            raise GateError(f"Invalid QUAST metric {name}: {fields[name]}") from error
        if not math.isfinite(number) or number < 0:
            raise GateError(f"Non-finite or negative QUAST metric {name}")
        numbers.append(number)
    return Metrics(*numbers)


@dataclass(frozen=True, slots=True)
class Telemetry:
    validated_k: int | None
    validated_k_minus_1: int | None
    rejected_ambiguous: int | None
    rejected_invalid: int | None
    active_nodes_after: int
    graph_conservation: str = "NOT_INSTRUMENTED_IN_ASSEMBLY;SEE_GRAPH_TESTS"


def telemetry(log: str) -> Telemetry:
    pattern = r"Junction stats: valid k joins = (\d+), valid k-1 joins = (\d+), ambiguous rejected = (\d+), invalid rejected = (\d+)"
    joins = re.findall(pattern, log)
    nodes = re.findall(r"Graph after exSPAnder: (\d+) active unitigs", log)
    if len(nodes) != 1:
        raise GateError("Missing or ambiguous post-exSPAnder node telemetry")
    if len(joins) > 1:
        raise GateError("Ambiguous junction telemetry")
    if joins:
        k, km1, amb, invalid = map(int, joins[0])
        return Telemetry(k, km1, amb, invalid, int(nodes[0]))
    return Telemetry(None, None, None, None, int(nodes[0]))


def acceptance(a: Metrics, b: Metrics, c: Metrics) -> tuple[str, ...]:
    reasons: list[str] = []
    for label, base in (("A", a), ("B", b)):
        for name in (
            "extensive_misassemblies",
            "local_misassemblies",
            "mismatches_per_100kbp",
            "indels_per_100kbp",
        ):
            if getattr(c, name) > getattr(base, name):
                reasons.append(
                    f"{name} increased vs {label}: {getattr(base, name)} -> {getattr(c, name)}"
                )
        if c.genome_fraction < base.genome_fraction:
            reasons.append(
                f"genome_fraction decreased vs {label}: {base.genome_fraction} -> {c.genome_fraction}"
            )
    return tuple(reasons)
