#!/usr/bin/env bash
#SBATCH --job-name=fastBact-F25-BM
#SBATCH --partition=qfree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=03:00:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/benchmark-results/f25_validation_20261030/logs/slurm-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/fastBact-phase-f25-correctness-worktree/benchmark-results/f25_validation_20261030/logs/slurm-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC RUST_LOG=info NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

P=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development
W_F25=$P/fastBact-phase-f25-correctness-worktree
R=$W_F25/benchmark-results/f25_validation_20261030
D=$P/04_ProkaDiff/data

BIN_CAND=$W_F25/target/release/fastBact
Q=/hpcfs/fhome/caizhh/.local/bin/quast.py

fail() { echo "BLOCKED: $*" >&2; exit 1; }

[[ -n "${SLURM_JOB_ID:-}" ]] || fail "SLURM_JOB_ID required"
[[ -x "$BIN_CAND" ]] || fail "Missing F25 binary: $BIN_CAND"
[[ -x "$Q" ]] || fail "Missing QUAST binary: $Q"

mkdir -p "$R/logs" "$R/assemblies" "$R/quast"
cd "$R"

echo "=== F25 BENCHMARK & CROSS-DATASET VALIDATION EXECUTION ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
echo "Slurm Job: $SLURM_JOB_ID, Node: $(hostname)"
echo "Candidate F25 SHA-256: $(sha256sum "$BIN_CAND" | awk '{print $1}')"

printf '%s\n' 'dataset	version	replicate	min_contig_bp	N50_bp	NGA50_bp	genome_fraction_pct	misassemblies	local_misassemblies	mismatches_per_100kbp	indels_per_100kbp	assembly_wall_s	assembly_peak_rss_kib	junction_valid_k	junction_valid_k_minus_1	junction_rejected_ambiguous	junction_rejected_invalid	active_unitigs_after_exspander' \
  > assembly_metrics.tsv.partial

extract() {
  awk -F '\t' -v OFS='\t' -v d="$1" -v v="$2" '
  function bad(s) {
    print "INVALID METRICS: " s > "/dev/stderr"
    failed=1
    exit 1
  }
  function num(s) { return s ~ /^[0-9]+([.][0-9]+)?$/ }
  BEGIN {
    n=split("N50|NGA50|Genome fraction (%)|# misassemblies|# local misassemblies|# mismatches per 100 kbp|# indels per 100 kbp",keys,"|")
    for(i=1;i<=n;i++) wanted[keys[i]]=1
  }
  FILENAME==ARGV[1] {
    if(NF!=2 || seen[$1]++) bad("malformed/duplicate QUAST row")
    if($1 in wanted) {
      if(!num($2)) bad($1)
      q[$1]=$2
    }
    next
  }
  FILENAME==ARGV[2] {
    if(NF!=3 || ++times!=1 || !num($1) || $1<=0 ||
       $2!~/^[0-9]+$/ || $2<=0 || $3!="0") bad("time/exit")
    wall=$1
    rss=$2
    next
  }
  /Junction stats:/ {
    line=$0
    sub(/^.*Junction stats: /,"",line)
    if(line!~/^valid k joins = [0-9]+, valid k-1 joins = [0-9]+, ambiguous rejected = [0-9]+, invalid rejected = [0-9]+[[:space:]]*$/)
      bad("junction syntax")
    if(split(line,j,/, /)!=4 || ++joins!=1)
      bad("junction cardinality")
    for(i=1;i<=4;i++) {
      sub(/^.*= /,"",j[i])
      sub(/[[:space:]]+$/,"",j[i])
    }
  }
  /Graph after exSPAnder:/ {
    line=$0
    sub(/^.*Graph after exSPAnder: /,"",line)
    if(line!~/^[0-9]+ active unitigs/) bad("node syntax")
    sub(/ .*/,"",line)
    nodes=line
    ++node_rows
  }
  END {
    if(failed) exit 1
    for(i=1;i<=n;i++)
      if(!(keys[i] in q)) bad("missing " keys[i])
    if(times!=1 || joins!=1 || node_rows!=1) bad("missing telemetry")
    if(q[keys[3]]>100 || q[keys[1]]<=0 || q[keys[2]]<=0) bad("range")
    printf "%s\t%s\t1\t200",d,v
    for(i=1;i<=n;i++) printf "\t%s",q[keys[i]]
    printf "\t%s\t%s",wall,rss
    for(i=1;i<=4;i++) printf "\t%s",j[i]
    printf "\t%s\n",nodes
  }' "$3" "$4" "$5"
}

run_dataset() {
  local ds=$1 v=$2 s=$3 ref=$4 bin=$BIN_CAND out q stem
  out=$R/assemblies/$ds/$v
  q=$R/quast/$ds/$v
  stem=$R/logs/$ds-$v
  rm -rf "$out" "$q"

  echo "--- Running $ds ($v) with sample $s ---"
  /usr/bin/time -f '%e\t%M\t%x' -o "$stem.assembly.time.tsv" \
    "$bin" \
    -1 "$D/prjna1088182/fastq/$s/"$s"_1.fastq.gz" \
    -2 "$D/prjna1088182/fastq/$s/"$s"_2.fastq.gz" \
    -o "$out" -t 8 --dump-stages "$out/stages" \
    > "$stem.assembly.stdout" 2> "$stem.assembly.stderr"

  [[ -s $out/contigs.fa ]] || fail "Missing $out/contigs.fa"

  cat "$stem.assembly.stdout" "$stem.assembly.stderr" |
    sed -E 's/\x1B\[[0-9;]*[mK]//g' > "$stem.assembly.combined.txt"

  echo "--- Running QUAST for $ds ($v) ---"
  /usr/bin/time -f '%e\t%M\t%x' -o "$stem.quast.time.tsv" \
    "$Q" -R "$ref" "$out/contigs.fa" -o "$q" \
    --threads 8 --min-contig 200 --silent \
    > "$stem.quast.stdout" 2> "$stem.quast.stderr"

  extract "$ds" "$v" "$q/report.tsv" "$stem.assembly.time.tsv" \
    "$stem.assembly.combined.txt" >> assembly_metrics.tsv.partial
  echo "--- Completed $ds ($v) ---"
}

# 1. Primary Acceptance Datasets
run_dataset R1-001 F25 SRR28349173 "$D/references/m6_identity/NC_000913.3.fasta"
run_dataset R1-003 F25 SRR28349164 "$D/references/m6_identity/CP009273.1.fasta"

# 2. Cross-Dataset Robustness Datasets (Independent biological samples)
run_dataset R1-001-Val F25 SRR28349174 "$D/references/m6_identity/NC_000913.3.fasta"
run_dataset R1-003-Val F25 SRR28349165 "$D/references/m6_identity/CP009273.1.fasta"

mv assembly_metrics.tsv.partial assembly_metrics.tsv
echo "=== F25 BENCHMARK & CROSS-DATASET VALIDATION COMPLETE ==="
date -u +"%Y-%m-%dT%H:%M:%SZ"
cat assembly_metrics.tsv
