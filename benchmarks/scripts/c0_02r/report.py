from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import platform
import sys

from .artifacts import Assembly, ComparedEvent, write_assemblies, write_events
from .execution import Command, Source
from .inputs import FileIdentity


@dataclass(frozen=True, slots=True)
class BenchmarkReport:
    status: str
    failures: list[str]
    sources: dict[str, Source]
    versions: dict[str, str]
    runs: list[Assembly]
    events: list[ComparedEvent]
    commands: list[Command]
    tools: list[FileIdentity]


def write_report(out: Path, report: BenchmarkReport) -> None:
    write_assemblies(out / "C0_02R_ASSEMBLY_AB.tsv", report.runs)
    write_events(out / "C0_02R_ALIGNMENT_EVENTS.tsv", report.events)
    write_events(
        out / "C0_02_JUNCTION_DELTAS.tsv",
        [e for e in report.events if e.comparison == "A_to_C"],
    )
    manifest = asdict(report)
    manifest.update(
        {
            "invocation": [sys.executable, *sys.argv],
            "node": platform.node(),
            "slurm": {
                key: os.environ.get(key)
                for key in (
                    "SLURM_JOB_ID",
                    "SLURMD_NODENAME",
                    "SLURM_JOB_PARTITION",
                    "SLURM_CPUS_PER_TASK",
                    "SLURM_MEM_PER_NODE",
                )
            },
            "environment": {
                key: os.environ.get(key)
                for key in ("CARGO_BUILD_JOBS", "RAYON_NUM_THREADS")
            },
            "coordinate_system": "0-based half-open; cs alleles follow target/aligned-query strand",
            "reference_correction_policy": "REFERENCE_COORDINATE_ONLY does not confirm correction; executed joins are not reference-confirmed events",
            "graph_inventory_scope": "Built-in stage dumps; after exSPAnder includes simplify; failure atomicity established separately by graph tests",
            "reproduction_scope": "REPRODUCIBLE_ON_ORIGINAL_HPC",
            "SELF_CONTAINED_EXTERNAL_REPRODUCTION": False,
        }
    )
    (out / "C0_02R_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    verdict = [
        f"# C0-02R benchmark acceptance: {report.status}",
        "",
        f"Completed assembly/QUAST pairs: {len(report.runs)}/18.",
        "",
        "Criteria: no increase in extensive/local misassemblies, mismatches or indels; no decrease in genome fraction versus either A or B. N50 is recorded and is not a sole acceptance criterion.",
        "",
        "Graph inventories use original built-in --dump-stages; after exSPAnder includes simplification. Counts do not establish path-level conservation; see deterministic G tests.",
        "",
        "Reference coordinates do not prove correction. Ambiguous repetitive mappings are explicitly retained. All coordinates are 0-based half-open; alleles follow alignment orientation.",
        "",
        "Reproduction: REPRODUCIBLE_ON_ORIGINAL_HPC. Inputs are not packaged for SELF_CONTAINED_EXTERNAL_REPRODUCTION.",
        "",
        *["- " + failure for failure in report.failures],
    ]
    (out / "C0_02R_FINAL_VERDICT.md").write_text("\n".join(verdict) + "\n")
