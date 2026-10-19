from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))
from c0_02r.alignment import AlignmentError, parse_paf, project_query
from c0_02r.execution import Executor, GateError
from c0_02r.metrics import QUAST_FIELDS, parse_quast
from c0_02r.structure import extract_events


def paf(cs: str, qlen: int, tlen: int) -> str:
    return f"q\t{qlen}\t0\t{qlen}\t+\tt\t{tlen}\t0\t{tlen}\t10\t{max(qlen, tlen)}\t60\ttp:A:P\tcs:Z:{cs}"


def test_deletion_consumes_target_coordinates_only() -> None:
    events = extract_events(parse_paf(paf(":5-a:5", 10, 11)))
    assert len(events) == 1
    e = events[0]
    assert (
        e.event_type,
        e.baseline_start,
        e.baseline_end,
        e.fixed_start,
        e.fixed_end,
    ) == ("DELETION", 5, 6, 5, 5)
    assert "CONFIRMED" not in e.evidence


def test_insertion_consumes_query_coordinates_only() -> None:
    e = extract_events(parse_paf(paf(":5+acc:5", 13, 10)))[0]
    assert (
        e.event_type,
        e.baseline_start,
        e.baseline_end,
        e.fixed_start,
        e.fixed_end,
    ) == ("INSERTION", 5, 5, 5, 8)


def test_reverse_substitution_uses_forward_query_coordinates() -> None:
    line = "q\t10\t0\t10\t-\tt\t110\t100\t110\t9\t10\t60\tcs:Z::2*ac:7"
    events = extract_events(parse_paf(line))
    e = next(e for e in events if e.event_type == "SUBSTITUTION")
    assert (e.baseline_start, e.fixed_start, e.fixed_end) == (102, 7, 8)
    assert events[0].event_type == "REVERSE_COMPLEMENT_ALIGNMENT"


def test_reverse_reference_projection_accounts_for_nonzero_query_start() -> None:
    a = parse_paf("q\t20\t5\t15\t-\tref\t200\t100\t110\t10\t10\t60\tcs:Z::10")[0]
    assert project_query(a, 5) == 109
    assert project_query(a, 14) == 100
    assert project_query(a, 15) is None


def test_contig_split_is_distinct_from_sequence_events() -> None:
    text = "q1\t5\t0\t5\t+\tt\t10\t0\t5\t5\t5\t60\tcs:Z::5\nq2\t5\t0\t5\t+\tt\t10\t5\t10\t5\t5\t60\tcs:Z::5"
    events = extract_events(parse_paf(text))
    assert [e.event_type for e in events] == ["CONTIG_SPLIT"]


def test_contig_merge_is_distinct_from_sequence_events() -> None:
    text = "q\t10\t0\t5\t+\tt1\t5\t0\t5\t5\t5\t60\tcs:Z::5\nq\t10\t5\t10\t+\tt2\t5\t0\t5\t5\t5\t60\tcs:Z::5"
    assert [e.event_type for e in extract_events(parse_paf(text))] == ["CONTIG_MERGE"]


def test_repetitive_secondary_alignment_does_not_duplicate_indel_events() -> None:
    primary = paf(":5-a:5", 10, 11)
    secondary = primary.replace("\tt\t", "\trepeat\t").replace("tp:A:P", "tp:A:S")
    events = extract_events(parse_paf(primary + "\n" + secondary))
    assert len([e for e in events if e.event_type == "DELETION"]) == 1
    assert all(e.evidence == "AMBIGUOUS_REPETITIVE_ALIGNMENT" for e in events)


@pytest.mark.parametrize("cs", [":5-a:4", ":5~gt10ag:5", ":5garbage:5"])
def test_malformed_or_inconsistent_alignment_fails_closed(cs: str) -> None:
    with pytest.raises(AlignmentError):
        parse_paf(paf(cs, 10, 11))


@pytest.mark.parametrize("bad", ["missing", "nan", "-", ""])
def test_quast_missing_or_unparseable_metrics_fail_closed(
    tmp_path: Path, bad: str
) -> None:
    report = tmp_path / "report.tsv"
    rows = [f"{field}\t1" for field in QUAST_FIELDS]
    if bad == "missing":
        rows.pop()
    else:
        rows[-1] = f"{QUAST_FIELDS[-1]}\t{bad}"
    report.write_text("\n".join(rows) + "\n")
    with pytest.raises(GateError):
        parse_quast(report)


def test_nonzero_command_is_recorded_and_rejected(tmp_path: Path) -> None:
    executor = Executor(tmp_path, [])
    with pytest.raises(GateError):
        executor.run((sys.executable, "-c", "raise SystemExit(17)"), "failure")
    assert executor.history[0].exit_code == 17
    assert (tmp_path / "failure.json").is_file()


def test_command_labels_with_periods_keep_distinct_receipts(tmp_path: Path) -> None:
    executor = Executor(tmp_path, [])
    first = executor.run((sys.executable, "-c", "print('A')"), "SYN-2.5M-A")
    second = executor.run((sys.executable, "-c", "print('B')"), "SYN-2.5M-B")
    assert first.stdout != second.stdout
    assert Path(first.stdout).read_text().strip() == "A"
    assert Path(second.stdout).read_text().strip() == "B"


def test_overlapping_split_is_reported_with_structural_ambiguity() -> None:
    text = "q1\t600\t0\t600\t+\tt\t1000\t0\t600\t600\t600\t60\ttp:A:P\tcs:Z::600\nq2\t600\t0\t600\t+\tt\t1000\t400\t1000\t600\t600\t60\ttp:A:P\tcs:Z::600"
    events = extract_events(parse_paf(text))
    assert len(events) == 1
    assert events[0].event_type == "CONTIG_SPLIT"
    assert "AMBIGUOUS" in events[0].evidence


def test_overlapping_merge_is_reported_with_structural_ambiguity() -> None:
    text = "q\t1000\t0\t600\t+\tt1\t600\t0\t600\t600\t600\t60\ttp:A:P\tcs:Z::600\nq\t1000\t400\t1000\t+\tt2\t600\t0\t600\t600\t600\t60\ttp:A:P\tcs:Z::600"
    events = extract_events(parse_paf(text))
    assert any(e.event_type == "CONTIG_MERGE" for e in events)
    assert any("AMBIGUOUS" in e.evidence for e in events)


def test_primary_inversion_tag_is_preserved() -> None:
    text = paf(":10", 10, 10).replace("\t+\t", "\t-\t").replace("tp:A:P", "tp:A:I")
    assert [e.event_type for e in extract_events(parse_paf(text))] == [
        "REVERSE_COMPLEMENT_ALIGNMENT"
    ]
    assert extract_events(parse_paf(text.replace("tp:A:I", "tp:A:i"))) == []


def test_unknown_alignment_type_fails_closed() -> None:
    with pytest.raises(AlignmentError):
        parse_paf(paf(":10", 10, 10).replace("tp:A:P", "tp:A:Q"))


def test_orchestrator_usage_loads_all_modules_without_computation() -> None:
    from c0_02r.benchmark import main

    assert main([]) == 2


def test_failure_report_keeps_fail_status_and_writes_isolated_artifacts(
    tmp_path: Path,
) -> None:
    import json
    from c0_02r.report import BenchmarkReport, write_report

    write_report(
        tmp_path, BenchmarkReport("FAIL", ["QUAST rejected"], {}, {}, [], [], [], [])
    )
    manifest = json.loads((tmp_path / "C0_02R_RUN_MANIFEST.json").read_text())
    assert manifest["status"] == "FAIL"
    assert manifest["failures"] == ["QUAST rejected"]
    assert (tmp_path / "C0_02_JUNCTION_DELTAS.tsv").read_text().count("\n") == 1


def test_missing_graph_stage_fails_closed(tmp_path: Path) -> None:
    from c0_02r.inputs import graph_stages

    path = tmp_path / "graph.json"
    path.write_text('{"stages": []}')
    with pytest.raises(GateError):
        graph_stages(path)


def test_n50_decrease_alone_does_not_fail_acceptance() -> None:
    from c0_02r.metrics import Metrics, acceptance

    a = Metrics(1000, 900, 99.9, 0, 0, 1, 1)
    c = Metrics(900, 900, 99.9, 0, 0, 1, 1)
    assert acceptance(a, a, c) == ()
