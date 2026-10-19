from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from c0_02r.execution import GateError
from c0_02r.verification import verify_manifest


@pytest.mark.parametrize(
    "text",
    [
        "{}",
        '{"status":"FAIL"}',
        '{"status":"PASS","runs":[],"commands":[],"sources":{}}',
    ],
)
def test_incomplete_or_failed_manifest_is_rejected(tmp_path: Path, text: str) -> None:
    path = tmp_path / "manifest.json"
    path.write_text(text)
    with pytest.raises(GateError):
        verify_manifest(path)


def manifest_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    from dataclasses import asdict
    import json
    from c0_02r import verification
    from c0_02r.execution import Source, digest
    from c0_02r.metrics import Metrics, QUAST_FIELDS

    source = Source(str(tmp_path), "c" * 40, "t" * 40, "s" * 64, "")

    def source_for(path: Path) -> Source:
        assert path == tmp_path
        return source

    monkeypatch.setattr(verification, "source_identity", source_for)
    data = tmp_path / "data.txt"
    data.write_text("fixture\n")
    file = {"path": str(data), "sha256": digest(data)}
    quast = tmp_path / "quast.tsv"
    quast.write_text("\n".join(f"{name}\t1" for name in QUAST_FIELDS) + "\n")
    runs = []
    for name in (
        "S1-COV-150X",
        "S1-COV-003",
        "SYN-2.5M",
        "S2-REAL-001",
        "R1-001",
        "R1-003",
    ):
        for variant in ("A", "B", "C"):
            runs.append(
                {
                    "dataset": {
                        "name": name,
                        "r1": file,
                        "r2": file,
                        "reference": file,
                    },
                    "variant": variant,
                    "binary": file,
                    "fasta": file,
                    "graph_telemetry": file,
                    "quast_report": {"path": str(quast), "sha256": digest(quast)},
                    "metrics": asdict(Metrics(1, 1, 1, 1, 1, 1, 1)),
                }
            )
    command = {
        "argv": ["fixture"],
        "cwd": str(tmp_path),
        "exit_code": 0,
        "stdout": str(data),
        "stderr": str(data),
        "stdout_sha256": digest(data),
        "stderr_sha256": digest(data),
    }
    manifest = {
        "status": "PASS",
        "runs": runs,
        "commands": [command],
        "sources": {name: asdict(source) for name in ("A", "B", "C")},
    }
    for filename in (
        "C0_02R_ASSEMBLY_AB.tsv",
        "C0_02R_ALIGNMENT_EVENTS.tsv",
        "C0_02_JUNCTION_DELTAS.tsv",
        "C0_02R_FINAL_VERDICT.md",
    ):
        (tmp_path / filename).write_text("fixture\n")
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    assert verify_manifest(path)[1] == 1
    return path


def test_recorded_metrics_must_match_hashed_quast_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import json

    path = manifest_fixture(tmp_path, monkeypatch)
    manifest = json.loads(path.read_text())
    manifest["runs"][0]["metrics"]["n50"] = 999
    path.write_text(json.dumps(manifest))
    with pytest.raises(GateError):
        verify_manifest(path)


def test_command_hash_cannot_override_other_artifact_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import json

    path = manifest_fixture(tmp_path, monkeypatch)
    manifest = json.loads(path.read_text())
    for run in manifest["runs"]:
        for field in ("r1", "r2", "reference"):
            run["dataset"][field]["sha256"] = "0" * 64
        for field in ("binary", "fasta", "graph_telemetry"):
            run[field]["sha256"] = "0" * 64
    path.write_text(json.dumps(manifest))
    with pytest.raises(GateError):
        verify_manifest(path)
