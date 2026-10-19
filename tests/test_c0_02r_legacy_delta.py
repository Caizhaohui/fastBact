from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))
from run_c0_02_benchmark import compare_fasta


def test_one_indel_does_not_create_coordinate_shift_substitutions(
    tmp_path: Path,
) -> None:
    rng = random.Random(20261019)
    dna = "".join(rng.choices("ACGT", k=12000))
    old, new = tmp_path / "old.fa", tmp_path / "new.fa"
    old.write_text(f">contig\n{dna[:6000]}A{dna[6000:]}\n")
    new.write_text(f">contig\n{dna}\n")
    events = compare_fasta(str(old), str(new))
    assert len(events) == 1
    assert events[0]["event_type"] == "DELETION"
    assert events[0]["event_len"] == 1


def test_reverse_complement_is_an_alignment_event(tmp_path: Path) -> None:
    rng = random.Random(2468)
    dna = "".join(rng.choices("ACGT", k=12000))
    reverse = dna.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    old, new = tmp_path / "old.fa", tmp_path / "new.fa"
    old.write_text(f">old\n{dna}\n")
    new.write_text(f">new\n{reverse}\n")
    events = compare_fasta(str(old), str(new))
    assert len(events) == 1
    assert events[0]["event_type"] == "REVERSE_COMPLEMENT_ALIGNMENT"
