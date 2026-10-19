from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Sequence
from .alignment import Event
from .execution import Command, Source
from .inputs import Dataset, FileIdentity, GraphStage
from .metrics import Metrics, Telemetry


@dataclass(frozen=True, slots=True)
class Assembly:
    dataset: Dataset
    variant: str
    source: Source
    binary: FileIdentity
    fasta: FileIdentity
    quast_report: FileIdentity
    graph_telemetry: FileIdentity
    stages: tuple[GraphStage, ...]
    metrics: Metrics
    junctions: Telemetry
    assembly_command: Command
    quast_command: Command


@dataclass(frozen=True, slots=True)
class ComparedEvent:
    dataset: str
    comparison: str
    event: Event


def write_events(path: Path, events: Sequence[ComparedEvent]) -> None:
    with path.open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset", "comparison", *[f.name for f in fields(Event)]])
        for row in events:
            writer.writerow([row.dataset, row.comparison, *asdict(row.event).values()])


def write_assemblies(path: Path, runs: Sequence[Assembly]) -> None:
    with path.open("w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "dataset",
                "kind",
                "variant",
                "commit",
                "source_tree_sha256",
                "dirty_status",
                "binary_sha256",
                "fasta_sha256",
                *[f.name for f in fields(Metrics)],
                "validated_k_joins",
                "validated_k_minus_1_joins",
                "rejected_ambiguous",
                "rejected_invalid",
                "nodes_before_exspander",
                "edges_before_exspander",
                "nodes_after_exspander",
                "edges_after_exspander",
                "graph_inventory_scope",
                "graph_telemetry_sha256",
                "peak_rss_kb",
                "wall_seconds",
                "assembly_exit_code",
                "quast_exit_code",
                "quast_report_sha256",
            ]
        )
        for run in runs:
            index = next(
                i for i, s in enumerate(run.stages) if s.name == "stage_05_exspander"
            )
            before, after = run.stages[index - 1 : index + 1]
            t = run.junctions
            writer.writerow(
                [
                    run.dataset.name,
                    run.dataset.kind,
                    run.variant,
                    run.source.commit,
                    run.source.source_tree_sha256,
                    run.source.dirty_status or "CLEAN",
                    run.binary.sha256,
                    run.fasta.sha256,
                    *asdict(run.metrics).values(),
                    t.validated_k,
                    t.validated_k_minus_1,
                    t.rejected_ambiguous,
                    t.rejected_invalid,
                    before.nodes,
                    before.edges,
                    after.nodes,
                    after.edges,
                    "BUILTIN_DUMP_STAGES;AFTER_INCLUDES_SIMPLIFY;SEE_ATOMICITY_TESTS",
                    run.graph_telemetry.sha256,
                    run.assembly_command.peak_rss_kb,
                    run.assembly_command.wall_seconds,
                    run.assembly_command.exit_code,
                    run.quast_command.exit_code,
                    run.quast_report.sha256,
                ]
            )
