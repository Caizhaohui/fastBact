from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from .execution import GateError, digest


@dataclass(frozen=True, slots=True)
class FileIdentity:
    path: str
    sha256: str


def identify(path: Path) -> FileIdentity:
    return FileIdentity(str(path.resolve(strict=True)), digest(path))


@dataclass(frozen=True, slots=True)
class Dataset:
    name: str
    kind: str
    r1: FileIdentity
    r2: FileIdentity
    reference: FileIdentity
    reproduction_scope: str = "REPRODUCIBLE_ON_ORIGINAL_HPC"
    self_contained_external_reproduction: bool = False


def load_datasets(config: Path, root: Path) -> tuple[Dataset, ...]:
    datasets: list[Dataset] = []
    with config.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if set(row) != {"dataset", "kind", "r1", "r2", "reference"}:
                raise GateError("Unexpected dataset TSV schema")
            datasets.append(
                Dataset(
                    row["dataset"],
                    row["kind"],
                    identify(root / row["r1"]),
                    identify(root / row["r2"]),
                    identify(root / row["reference"]),
                )
            )
    expected = {
        "S1-COV-150X",
        "S1-COV-003",
        "SYN-2.5M",
        "S2-REAL-001",
        "R1-001",
        "R1-003",
    }
    if len(datasets) != 6 or {d.name for d in datasets} != expected:
        raise GateError("Benchmark requires exactly the six requested datasets")
    return tuple(datasets)


@dataclass(frozen=True, slots=True)
class GraphStage:
    name: str
    nodes: int
    edges: int
    bases: int


def graph_stages(path: Path) -> tuple[GraphStage, ...]:
    match json.loads(path.read_text()):
        case {"stages": list(stages)}:
            parsed: list[GraphStage] = []
            for stage in stages:
                match stage:
                    case {
                        "stage_name": str(name),
                        "active_unitigs": int(nodes),
                        "active_edges": int(edges),
                        "total_graph_bp": int(bases),
                    }:
                        if min(nodes, edges, bases) < 0:
                            raise GateError("Negative graph inventory")
                        parsed.append(GraphStage(name, nodes, edges, bases))
                    case _:
                        raise GateError("Malformed graph telemetry stage")
            if (
                len(parsed) < 3
                or sum(s.name == "stage_05_exspander" for s in parsed) != 1
            ):
                raise GateError("Missing or duplicated exSPAnder stage telemetry")
            return tuple(parsed)
        case _:
            raise GateError("Missing graph telemetry stages")
