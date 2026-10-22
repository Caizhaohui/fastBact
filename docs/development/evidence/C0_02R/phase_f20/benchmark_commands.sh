#!/usr/bin/env bash
#SBATCH --job-name=fastBact-F20
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=02:30:00
#SBATCH --output=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmark-results/C0_02R_f20_selective_20261022/logs/slurm-%j.out
#SBATCH --error=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmark-results/C0_02R_f20_selective_20261022/logs/slurm-%j.err

set -Eeuo pipefail
umask 027
export LC_ALL=C TZ=UTC RUST_LOG=info NO_COLOR=1
export RAYON_NUM_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1

P=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development
REPO=$P/01_fastBact
W=$P/fastBact-phase-f19-codex-worktree
R=$REPO/benchmark-results/C0_02R_f20_selective_20261022
E=$REPO/docs/development/evidence/C0_02R/phase_f20
D=$P/04_ProkaDiff/data

A=$REPO/benchmark-results/C0_02R_pe_fix_20261019/targets/B/release/fastBact
B=$W/target/release/fastBact
Q=$(readlink -f /hpcfs/fhome/caizhh/.local/bin/quast.py)

fail() { echo "BLOCKED: $*" >&2; exit 1; }

[[ -n "$SLURM_JOB_ID" && $SLURM_JOB_PARTITION == qfree ]] ||
  fail "Slurm/qfree required"
[[ $SLURM_CPUS_PER_TASK == 8 && $SLURM_MEM_PER_NODE == 12288 ]] ||
  fail "8 CPUs/12G required"
[[ $(scontrol show hostnames "$SLURM_JOB_NODELIST") == vnode1 ]] ||
  fail "vnode1 required"

for t in sha256sum gzip jq awk srun git sed; do
  command -v "$t" >/dev/null || fail "Missing $t"
done
[[ -x /usr/bin/time && -x $A && -x $B && -x $Q ]] ||
  fail "Missing executable"

mkdir -p "$R/logs" "$E"
mkdir "$R/.execution-lock" || fail "Results root already claimed"
STATUS=BLOCKED

finish() {
  rc=$?
  trap - EXIT
  printf '%s\t%s\t%s\n' "$STATUS" "$rc" "$SLURM_JOB_ID" > "$R/status.tsv"
  cp "$R/status.tsv" "$E/run-$SLURM_JOB_ID-status.tsv" || rc=1
  exit "$rc"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

for x in assemblies quast assembly_metrics.tsv assembly_metrics.tsv.partial; do
  [[ ! -e $R/$x ]] || fail "Existing $x"
done

cd "$R"
exec 9>"logs/commands-$SLURM_JOB_ID.log"
export BASH_XTRACEFD=9
set -x

scontrol show job -dd "$SLURM_JOB_ID" > logs/allocation.txt
uname -a > logs/host.txt
lscpu > logs/lscpu.txt
/usr/bin/time --version > logs/time-version.txt
"$Q" --version > logs/quast-version.txt
[[ $(cat logs/quast-version.txt) == 'QUAST v5.2.0' ]] ||
  fail "QUAST version drift"
sha256sum "$Q" > logs/quast-launcher.sha256
"$A" --help > logs/A-help.txt
"$B" --help > logs/B-help.txt

identity() {
  [[ $(git -C "$W" rev-parse HEAD) == "0606be423d187771f5925a83da2cc377c40c2a3d" ]] || fail "B HEAD"
  [[ $(git -C "$W" rev-parse HEAD:src) == "30967bf4940844aeb6b744732bcad6d2a72daab4" ]] || fail "B src tree"
  [[ -z $(git -C "$W" status --porcelain --untracked-files=all -- src Cargo.toml Cargo.lock) ]] || fail "Dirty production/build inputs"
}
identity

{
  printf '%s  %s\n' \
    2b184db0f108cb508ca7e591240eddb551bea8ebad90ac5eb8ee360bac78ea38 "$A"
  printf '%s  %s\n' \
    4fbafeeeffe94abeef42efaa5710f75dbd0af8b234edd2bffba43cc4f50198f0 "$B"

  for lock in "$W/Cargo.lock" \
    "$REPO/benchmark-results/C0_02R_pe_fix_20261019/sources/B/Cargo.lock"; do
    printf '%s  %s\n' \
      b4f683cbab9eabc7a484ae459641372082b486dccb26b90b14fa840bf4d40295 "$lock"
  done

  printf '%s  %s\n' \
    e25760e43d1126b9b56d4ee9de79a42cfd28b1de1dd54224bcc24319bbbe5ebf \
    "$D/prjna1088182/fastq/SRR28349173/SRR28349173_1.fastq.gz"
  printf '%s  %s\n' \
    ccf7c2627919cf2df6eb03c06c322d8dea241ad0ce9bdc34e7c6c13f4b983701 \
    "$D/prjna1088182/fastq/SRR28349173/SRR28349173_2.fastq.gz"
  printf '%s  %s\n' \
    dd43192faf592582c30dfcf71938d6450f6525b54ca254276dbd8ea16672011b \
    "$D/prjna1088182/fastq/SRR28349164/SRR28349164_1.fastq.gz"
  printf '%s  %s\n' \
    94fd5d094315c46725840a5ce8abbff4fc4842df748930b5b8f4bff9497fc3d1 \
    "$D/prjna1088182/fastq/SRR28349164/SRR28349164_2.fastq.gz"
  printf '%s  %s\n' \
    6b195feda4c66140f6762742eb8b30c2652f02b45878b174f5b00ef85ecc95d7 \
    "$D/references/m6_identity/NC_000913.3.fasta"
  printf '%s  %s\n' \
    37eade01311f1f81659fad53bf907a93139075a75d115fd3aec2ef09beac2e90 \
    "$D/references/m6_identity/CP009273.1.fasta"
} > logs/expected.sha256

sha256sum -c logs/expected.sha256 > logs/checksums-before.txt

for s in SRR28349173 SRR28349164; do
  gzip -t \
    "$D/prjna1088182/fastq/$s/"$s"_1.fastq.gz" \
    "$D/prjna1088182/fastq/$s/"$s"_2.fastq.gz"
done

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
    for(i=1;i<=2;i++)
      if(q[keys[i]]!=int(q[keys[i]])) bad("noninteger length")
    for(i=4;i<=5;i++)
      if(q[keys[i]]!=int(q[keys[i]])) bad("noninteger count")
    if(rss>12582912) bad("RSS exceeds allocation")
    printf "%s\t%s\t1\t200",d,v
    for(i=1;i<=n;i++) printf "\t%s",q[keys[i]]
    printf "\t%s\t%s",wall,rss
    for(i=1;i<=4;i++) printf "\t%s",j[i]
    printf "\t%s\n",nodes
  }' "$3" "$4" "$5"
}

run_one() {
  local ds=$1 v=$2 s=$3 ref=$4 bin out q stem
  bin=$A
  [[ $v == A ]] || bin=$B
  out=$R/assemblies/$ds/$v
  q=$R/quast/$ds/$v
  stem=$R/logs/$ds-$v

  [[ ! -e $out && ! -e $q ]] || fail "Run output exists"
  uptime > "$stem.load-before.txt"

  srun --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
    /usr/bin/time -f '%e\t%M\t%x' -o "$stem.assembly.time.tsv" \
    "$bin" \
    -1 "$D/prjna1088182/fastq/$s/"$s"_1.fastq.gz" \
    -2 "$D/prjna1088182/fastq/$s/"$s"_2.fastq.gz" \
    -o "$out" -t 8 --dump-stages "$out/stages" \
    > "$stem.assembly.stdout" 2> "$stem.assembly.stderr"

  [[ -s $out/contigs.fa && -s $out/stages/stage_telemetry.json ]] ||
    fail "Missing assembly artifacts"

  awk '
    /^>/ {
      if(headers && bases==0) exit 1
      headers++
      bases=0
      next
    }
    {
      if(!headers || $0!~/^[ACGTNacgtn]+$/) exit 1
      bases+=length($0)
      total+=length($0)
    }
    END { if(!headers || !bases || !total) exit 1 }
  ' "$out/contigs.fa" || fail "Invalid FASTA"

  jq -e '.subsampling_enabled == false and
    (.input_read_pairs > 0) and (.retained_read_pairs > 0) and
    (.retained_read_pairs <= .input_read_pairs) and
    (.stages | length >= 3) and
    (all(.stages[];
      .active_unitigs >= 0 and
      .active_edges >= 0 and
      .total_graph_bp >= 0)) and
    ([.stages[] | select(.stage_name == "stage_05_exspander")] |
      length == 1) and
    ([.stages[] | select(.stage_name == "stage_06_final_polished")] |
      length == 1)' \
    "$out/stages/stage_telemetry.json" > "$stem.stage-check.txt"

  cat "$stem.assembly.stdout" "$stem.assembly.stderr" |
    sed -E 's/\x1B\[[0-9;]*[mK]//g' > "$stem.assembly.combined.txt"

  srun --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
    /usr/bin/time -f '%e\t%M\t%x' -o "$stem.quast.time.tsv" \
    "$Q" -R "$ref" "$out/contigs.fa" -o "$q" \
    --threads 8 --min-contig 200 --silent \
    > "$stem.quast.stdout" 2> "$stem.quast.stderr"

  extract "$ds" "$v" "$q/report.tsv" "$stem.assembly.time.tsv" \
    "$stem.assembly.combined.txt" >> assembly_metrics.tsv.partial

  uptime > "$stem.load-after.txt"
}

run_one R1-001 A SRR28349173 "$D/references/m6_identity/NC_000913.3.fasta"
run_one R1-001 B SRR28349173 "$D/references/m6_identity/NC_000913.3.fasta"
run_one R1-003 B SRR28349164 "$D/references/m6_identity/CP009273.1.fasta"
run_one R1-003 A SRR28349164 "$D/references/m6_identity/CP009273.1.fasta"

sha256sum -c logs/expected.sha256 > logs/checksums-after.txt
sha256sum -c logs/quast-launcher.sha256 > logs/quast-checksum-after.txt
identity

awk -F '\t' -v OFS='\t' '
  NR==1 { next }
  {
    if(NF!=18 ||
       ($1!="R1-001" && $1!="R1-003") ||
       ($2!="A" && $2!="B") ||
       seen[$1,$2]++) {
      invalid=1
      exit 1
    }
    for(i=5;i<=11;i++) m[$1,$2,i]=$i
    n++
  }
  END {
    if(invalid || n!=4) exit 1
    print "dataset","metric_column","A","B","result"
    split("R1-001 R1-003",d," ")
    for(k=1;k<=2;k++) {
      if(!seen[d[k],"A"] || !seen[d[k],"B"]) exit 1
      for(i=5;i<=11;i++) {
        a=m[d[k],"A",i]+0
        b=m[d[k],"B",i]+0
        bad=(i<=7 ? b<a : b>a)
        print d[k],i,a,b,(bad ? "REGRESSION" : "PASS")
      }
    }
  }' assembly_metrics.tsv.partial > acceptance.tsv || {
    STATUS=REJECTED_OR_INVALID
    fail "Acceptance evaluation script failed"
  }

mv assembly_metrics.tsv.partial assembly_metrics.tsv
cp assembly_metrics.tsv acceptance.tsv "$E/"

find assemblies quast logs -type f \
  ! -name "slurm-*" ! -name "commands-*" -print0 |
  sort -z | xargs -0 sha256sum > artifact_checksums.sha256
sha256sum assembly_metrics.tsv acceptance.tsv >> artifact_checksums.sha256
cp artifact_checksums.sha256 "$E/"

if grep -q "REGRESSION" acceptance.tsv; then
  STATUS=REGRESSION_DETECTED
  echo "BENCHMARK FINISHED: REGRESSION DETECTED (see acceptance.tsv)" >&2
else
  STATUS=PASS_SELECTIVE
  echo "BENCHMARK FINISHED: PASS_SELECTIVE (all safety & contiguity criteria met)" >&2
fi
