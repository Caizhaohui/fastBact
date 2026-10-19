from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
import sys

from .alignment import AlignmentError
from .artifacts import Assembly, ComparedEvent
from .comparison import align, annotate_reference, minimap2, unaligned_events
from .execution import Executor, GateError, digest, require_compute, source_identity
from .inputs import FileIdentity, graph_stages, identify, load_datasets
from .report import BenchmarkReport, write_report
from .metrics import acceptance, parse_quast, telemetry
from .structure import extract_events


def main(arguments: list[str]) -> int:
    if len(arguments) != 3:
        print(
            "Usage: run_c0_02_benchmark.py SOURCES_ROOT DATA_ROOT NEW_OUTPUT_DIR",
            file=sys.stderr,
        )
        return 2
    require_compute()
    sources_root, data_root, out = (Path(arg).resolve() for arg in arguments)
    sources = {name: source_identity(sources_root / name) for name in ("A", "B", "C")}
    for source in sources.values():
        if source.dirty_status:
            raise GateError(
                f"Refusing dirty benchmark source: {source.path}: {source.dirty_status}"
            )
    if (
        sources["A"].commit != "3950b2ab2584b24d39172b25bbbcb565232891b7"
        or sources["B"].commit != "b9727c3ab372f7b995e439d92b72c7a4bfd051cf"
    ):
        raise GateError("Baseline A/B commit identity mismatch")
    out.mkdir(parents=True, exist_ok=False)
    executor = Executor(out / "commands", [])
    runs: list[Assembly] = []
    events: list[ComparedEvent] = []
    failures: list[str] = []
    versions: dict[str, str] = {
        "python": sys.version,
        "python_executable": sys.executable,
    }
    status = "FAIL"
    tools: list[FileIdentity] = []
    quast = os.environ.get("QUAST", str(Path.home() / ".local/bin/quast.py"))
    os.environ["CARGO_BUILD_JOBS"] = "4"
    os.environ["RAYON_NUM_THREADS"] = "8"
    try:
        tools.extend((identify(Path(quast)), identify(Path(minimap2()))))
        for label, command in (
            ("rust", ("rustc", "--version")),
            ("cargo", ("cargo", "--version")),
            ("quast", (quast, "--version")),
            ("minimap2", (minimap2(), "--version")),
            (
                "slurm_allocation",
                ("scontrol", "show", "job", "-o", os.environ["SLURM_JOB_ID"]),
            ),
        ):
            receipt = executor.run(command, "version-" + label)
            versions[label] = (
                Path(receipt.stdout).read_text().strip()
                + Path(receipt.stderr).read_text().strip()
            )
        c_root = Path(sources["C"].path)
        os.environ["CARGO_TARGET_DIR"] = str(sources_root.parent / "targets" / "C")
        executor.run(
            ("cargo", "test", "--release", "--locked", "--offline", "--all-targets"),
            "rust-regression",
            c_root,
        )
        executor.run(
            (
                sys.executable,
                "-m",
                "pytest",
                "tests/test_c0_02r_alignment.py",
                "tests/test_c0_02r_legacy_delta.py",
                "tests/test_benchmark_parser.py",
                "-v",
            ),
            "python-regression",
            c_root,
        )
        binaries = {}
        for variant, source in sources.items():
            target = sources_root.parent / "targets" / variant
            os.environ["CARGO_TARGET_DIR"] = str(target)
            executor.run(
                (
                    "cargo",
                    "build",
                    "--release",
                    "--locked",
                    "--offline",
                    "--bin",
                    "fastBact",
                ),
                "build-" + variant,
                Path(source.path),
            )
            binaries[variant] = identify(target / "release/fastBact")
        datasets = load_datasets(
            c_root / "benchmarks/configs/c0_02r_datasets.tsv", data_root
        )
        for dataset in datasets:
            print(f"RUNNING {dataset.name}", flush=True)
            current: dict[str, Assembly] = {}
            for variant, source in sources.items():
                directory = out / "assemblies" / dataset.name / variant
                assembly = executor.run(
                    (
                        binaries[variant].path,
                        "-1",
                        dataset.r1.path,
                        "-2",
                        dataset.r2.path,
                        "-o",
                        str(directory),
                        "-t",
                        "8",
                        "--dump-stages",
                        str(directory / "stages"),
                    ),
                    f"{dataset.name}-{variant}-assembly",
                )
                fasta = identify(directory / "contigs.fa")
                quast_run = executor.run(
                    (
                        quast,
                        "-R",
                        dataset.reference.path,
                        fasta.path,
                        "-o",
                        str(directory / "quast"),
                        "--threads",
                        "8",
                        "--min-contig",
                        "200",
                        "--silent",
                    ),
                    f"{dataset.name}-{variant}-quast",
                )
                report = identify(directory / "quast/report.tsv")
                graph = identify(directory / "stages/stage_telemetry.json")
                log = (
                    Path(assembly.stdout).read_text()
                    + Path(assembly.stderr).read_text()
                )
                joins = telemetry(log)
                if variant != "A" and joins.validated_k is None:
                    raise GateError(
                        f"Missing junction telemetry for {variant}/{dataset.name}"
                    )
                run = Assembly(
                    dataset,
                    variant,
                    source,
                    binaries[variant],
                    fasta,
                    report,
                    graph,
                    graph_stages(Path(graph.path)),
                    parse_quast(Path(report.path)),
                    joins,
                    assembly,
                    quast_run,
                )
                runs.append(run)
                current[variant] = run
            refs = {
                variant: align(
                    executor,
                    (Path(dataset.reference.path), Path(run.fasta.path)),
                    f"{dataset.name}-{variant}-reference",
                )
                for variant, run in current.items()
            }
            for variant, mappings in refs.items():
                for event in extract_events(mappings):
                    event = replace(
                        event,
                        ref_coords=f"{event.baseline_contig}:{event.baseline_start}-{event.baseline_end}",
                        evidence="REFERENCE_" + event.evidence,
                    )
                    events.append(
                        ComparedEvent(dataset.name, variant + "_vs_REFERENCE", event)
                    )
            for baseline, fixed in (("A", "B"), ("A", "C"), ("B", "C")):
                pair = (
                    Path(current[baseline].fasta.path),
                    Path(current[fixed].fasta.path),
                )
                mappings = align(
                    executor, pair, f"{dataset.name}-{baseline}-vs-{fixed}"
                )
                for event in extract_events(mappings) + unaligned_events(
                    pair, mappings
                ):
                    events.append(
                        ComparedEvent(
                            dataset.name,
                            baseline + "_to_" + fixed,
                            annotate_reference(event, refs[fixed]),
                        )
                    )
            failures.extend(
                dataset.name + ": " + reason
                for reason in acceptance(
                    current["A"].metrics, current["B"].metrics, current["C"].metrics
                )
            )
            if (
                dataset.name in {"S1-COV-150X", "S1-COV-003"}
                and len({run.fasta.sha256 for run in current.values()}) != 1
            ):
                failures.append(dataset.name + ": invariant FASTA identity changed")
            for identity in (dataset.r1, dataset.r2, dataset.reference):
                if digest(Path(identity.path)) != identity.sha256:
                    raise GateError(f"Input changed during benchmark: {identity.path}")
            print(
                f"COMPLETED {dataset.name}: cumulative acceptance failures={len(failures)}",
                flush=True,
            )
        for source in sources.values():
            if source_identity(Path(source.path)) != source:
                raise GateError(f"Source changed during benchmark: {source.path}")
        status = "PASS" if len(runs) == 18 and not failures else "FAIL"
    except (GateError, AlignmentError, OSError, ValueError) as error:
        failures.append(str(error))
        print(f"FAIL: {error}", file=sys.stderr, flush=True)
    finally:
        write_report(
            out,
            BenchmarkReport(
                status,
                failures,
                sources,
                versions,
                runs,
                events,
                executor.history,
                tools,
            ),
        )
    return 0 if status == "PASS" else 1
