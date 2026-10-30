#!/usr/bin/env bash
#SBATCH --job-name=fastBact-PostInt
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/slurm-postint-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/slurm-postint-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

cd /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact

echo "=== FASTBACT F25 POST-INTEGRATION VERIFICATION ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
echo "Slurm Job: $SLURM_JOB_ID, Node: $(hostname)"
echo "Git HEAD: $(git rev-parse HEAD)"

# 1. Full Release Test Suite
echo "--- Running Full Cargo Release Test Suite ---"
cargo test --release -- --nocapture

# 2. Release Binary Compilation
echo "--- Building Release Binary ---"
cargo build --release
BIN=target/release/fastBact
echo "Binary SHA-256: $(sha256sum "$BIN" | awk '{print $1}')"

# 3. CLI Smoke Test
echo "--- CLI Smoke Test ---"
"$BIN" --help | grep -q "fastBact"
"$BIN" --version

# 4. Synthetic Assembly Smoke Test
echo "--- Synthetic Assembly Smoke Test ---"
TMPDIR=$(mktemp -d /tmp/fastbact-postint-XXXXXX)
python3 -c "
import random
random.seed(42)
genome = ''.join(random.choices('ACGT', k=10000))
with open('$TMPDIR/ref.fa', 'w') as f:
    f.write('>ref\n' + genome + '\n')
reads1 = []
reads2 = []
for i in range(5000):
    pos = random.randint(0, len(genome) - 300)
    frag = genome[pos:pos+300]
    r1 = frag[:150]
    rc = {'A':'T','C':'G','G':'C','T':'A'}
    r2 = ''.join(rc.get(b, 'N') for b in reversed(frag[-150:]))
    reads1.append(f'@read_{i}/1\n{r1}\n+\n' + 'I'*150)
    reads2.append(f'@read_{i}/2\n{r2}\n+\n' + 'I'*150)
import gzip
with gzip.open('$TMPDIR/r1.fq.gz', 'wt') as f:
    f.write('\n'.join(reads1) + '\n')
with gzip.open('$TMPDIR/r2.fq.gz', 'wt') as f:
    f.write('\n'.join(reads2) + '\n')
"
"$BIN" -1 "$TMPDIR/r1.fq.gz" -2 "$TMPDIR/r2.fq.gz" -o "$TMPDIR/out" -t 4
[[ -s "$TMPDIR/out/contigs.fa" ]]
CONTIG_LEN=$(grep -v '^>' "$TMPDIR/out/contigs.fa" | tr -d '\n' | wc -c)
echo "Synthetic assembled bases: $CONTIG_LEN / 10000"
rm -rf "$TMPDIR"

echo "=== ALL POST-INTEGRATION CHECKS PASSED ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
