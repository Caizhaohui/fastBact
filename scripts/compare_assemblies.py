import sys

def analyze(fasta_path, name):
    lens = []
    seq = ""
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    lens.append(len(seq))
                    seq = ""
            else:
                seq += line
        if seq:
            lens.append(len(seq))
    lens.sort(reverse=True)
    total_len = sum(lens)
    half = total_len / 2
    acc = 0
    n50 = 0
    for l in lens:
        acc += l
        if acc >= half:
            n50 = l
            break
    max_len = lens[0] if lens else 0
    num_contigs = len(lens)
    plasmid_cand = [l for l in lens if abs(l - 20000) <= 200]
    print(f"=== {name} ===")
    print(f"Total contigs: {num_contigs}")
    print(f"Total length: {total_len} bp")
    print(f"Longest contig: {max_len} bp")
    print(f"N50: {n50} bp")
    print(f"Top 5 lengths: {lens[:5]}")
    print(f"Plasmid 20kb candidates: {plasmid_cand}")

analyze("benchmark/out_fastbact/contigs.fa", "fastBact")
analyze("benchmark/out_shovill/contigs.fa", "Shovill")
analyze("benchmark/out_spades/contigs.fasta", "SPAdes")
