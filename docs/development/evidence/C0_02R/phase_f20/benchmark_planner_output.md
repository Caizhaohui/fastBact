I could not save `docs/development/evidence/C0_02R/phase_f20/benchmark_planner_output.md`: this session enforces read-only filesystem access. The complete plan follows. No benchmark jobs were submitted.

Both binary hashes, both Cargo.lock hashes, all six input hashes, candidate HEAD, and candidate `src/` tree match the supplied identities. QUAST is version 5.2.0. Slurm controller access is blocked in this session.

---

**fastBact Phase F20 selective real-data benchmark plan**

**Scope and frozen identities**

Evaluate only R1-001 and R1-003, using the supplied paired FASTQs and corresponding references.

| F20 label | Source identity | Binary SHA-256 |
|---|---|---|
| A | `b9727c3ab372f7b995e439d92b72c7a4bfd051cf` | `2b184db0f108cb508ca7e591240eddb551bea8ebad90ac5eb8ee360bac78ea38` |
| B | `0606be423d187771f5925a83da2cc377c40c2a3d` | `4fbafeeeffe94abeef42efaa5710f75dbd0af8b234edd2bffba43cc4f50198f0` |

Candidate code commit: `a9930471a73612159a83694d3d661f060042efc8`.

Candidate `src/` tree: `30967bf4940844aeb6b744732bcad6d2a72daab4`.

Both supplied build records identify `rustc 1.89.0 (29483883e 2025-08-04)`. The frozen binary hashes are the execution identities; checking the current compiler would not independently establish their build provenance. Do not rebuild either binary.

F20 A corresponds to historical B in the earlier three-version harness. Do not invoke that harness unchanged: it rebuilds binaries, expects three variants and six datasets, and uses different labels.

**Execution design**

Use one Slurm allocation on `qfree/vnode1`, one task, eight CPUs, 12 GiB memory, six-hour job limit. Run assemblies and QUAST serially:

1. R1-001 A.
2. R1-001 B.
3. R1-003 B.
4. R1-003 A.

Reversing order for the second dataset partially balances execution order. One measurement per cell supports descriptive wall-time/RSS comparisons only; it does not support statistical performance claims.

Use the exact requested assembly options. Preserve defaults: automatic k-mers, polishing enabled, minimum output contig length 500 bp, minimum coverage 2, and automatic subsampling disabled. The default depth value alone does not enable subsampling. Do not add `--min-len 200` to match QUAST.

Coverage labels are approximate metadata. Record actual ingest counts and base counts from stage telemetry; do not use nominal coverage to change parameters.

**Exact `benchmark_commands.sh`**

Save this script under the evidence root. Its Bash syntax was checked. The extraction function was exercised against an archived real report/log and rejected an invalid-metric fixture. No assembly was executed during planning.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=fastBact-F20
#SBATCH --partition=qfree
#SBATCH --nodelist=vnode1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=06:00:00
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
  [[ $(git -C "$W" rev-parse HEAD) ==
     0606be423d187771f5925a83da2cc377c40c2a3d ]] || fail "B HEAD"
  [[ $(git -C "$W" rev-parse HEAD:src) ==
     30967bf4940844aeb6b744732bcad6d2a72daab4 ]] || fail "B src tree"
  [[ -z $(git -C "$W" status --porcelain --untracked-files=all \
    -- src Cargo.toml Cargo.lock) ]] || fail "Dirty production/build inputs"
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
        if(bad) rejected=1
      }
    }
    if(rejected) exit 2
  }' assembly_metrics.tsv.partial > acceptance.tsv || {
    STATUS=REJECTED_OR_INVALID
    fail "Acceptance failed; retain partial metrics and evidence"
  }

mv assembly_metrics.tsv.partial assembly_metrics.tsv
cp assembly_metrics.tsv acceptance.tsv "$E/"

find assemblies quast logs -type f \
  ! -name "slurm-*" ! -name "commands-*" -print0 |
  sort -z | xargs -0 sha256sum > artifact_checksums.sha256
sha256sum assembly_metrics.tsv acceptance.tsv >> artifact_checksums.sha256
cp artifact_checksums.sha256 "$E/"

STATUS=PASS_SELECTIVE
```

**Submission and environment capture**

Before submission, save the plan and script, create `logs/`, and verify that assembly and QUAST output directories do not already exist. The existing empty results root is permitted.

```bash
P=/hpcfs/fhome/caizhh/Desktop/03_Tool_Development
R=$P/01_fastBact/benchmark-results/C0_02R_f20_selective_20261022
E=$P/01_fastBact/docs/development/evidence/C0_02R/phase_f20

mkdir -p "$R/logs"
bash -n "$E/benchmark_commands.sh"
sha256sum "$E/benchmark_commands.sh" > "$E/benchmark_commands.sha256"
sbatch --parsable "$E/benchmark_commands.sh" > "$E/slurm_submission.txt"
```

Before the first measured assembly, the executor must also archive the QUAST installation identity: resolved launcher, interpreter/shebang, package installation location, installed package/version listing, and checksums of QUAST code and bundled alignment executables. Preserve this environment for the entire run. The script’s launcher hash alone does not fingerprint the complete QUAST installation.

Record available disk space and confirm adequate capacity before submission. Never delete existing evidence to make space.

**Metric extraction contract**

`assembly_metrics.tsv` is UTF-8, tab-delimited, with one header and exactly four data rows. Column order is defined by the script.

| Column | Source and meaning |
|---|---|
| `N50_bp` | QUAST `N50`, bp |
| `NGA50_bp` | QUAST `NGA50`, bp; never substitute N50 |
| `genome_fraction_pct` | QUAST `Genome fraction (%)`, percentage points |
| `misassemblies` | QUAST `# misassemblies` |
| `local_misassemblies` | QUAST `# local misassemblies` |
| `mismatches_per_100kbp` | QUAST `# mismatches per 100 kbp` |
| `indels_per_100kbp` | QUAST `# indels per 100 kbp` |
| `assembly_wall_s` | GNU time `%e`, seconds, assembly command only |
| `assembly_peak_rss_kib` | GNU time `%M`, KiB, assembly process |
| Four junction counters | Exact `Junction stats:` record |
| `active_unitigs_after_exspander` | Exact post-exSPAnder log record |

Timing excludes preflight checks and QUAST. QUAST timing is preserved separately. Slurm accounting is a separate resource receipt and must not replace GNU time values silently.

Both binaries must emit exactly one junction record and one post-exSPAnder node record. Parse stdout and stderr together: archived baseline telemetry appears in stdout. Missing telemetry is an error, never zero or `NA`.

Retain raw `report.tsv`, QUAST logs and alignment reports, both assembly streams, stage FASTAs, and `stage_telemetry.json`. Preserve ingest counts, base counts, subsampling policy, and graph inventories. Graph inventories and junction totals do not prove graph conservation; retain the qualification `NOT_INSTRUMENTED_IN_ASSEMBLY;SEE_GRAPH_TESTS`.

Missing, duplicate, malformed, negative, nonfinite, or nonnumeric required metrics invalidate the comparison. Undefined NGA50 is a blocked quality assessment, never a passing result.

**Acceptance boundaries**

Apply every comparison independently to each dataset, using the displayed QUAST values without an added tolerance.

| Gate | Required B versus A |
|---|---|
| Safety: genome fraction | B ≥ A |
| Safety: misassemblies | B ≤ A |
| Safety: local misassemblies | B ≤ A |
| Safety: mismatches/100 kbp | B ≤ A |
| Safety: indels/100 kbp | B ≤ A |
| Contiguity: N50 | B ≥ A |
| Contiguity: NGA50 | B ≥ A |

The safety boundaries preserve the existing benchmark policy. The N50/NGA50 boundaries are an additional conservative F20 contiguity screen. Report safety and contiguity outcomes separately; the script requires both for `PASS_SELECTIVE`.

Do not average away a regression across datasets. Improved N50 cannot compensate for increased errors or reduced genome fraction.

Report wall-time and RSS differences as absolute values and `100 × (B/A − 1)`. There is no statistical speed or memory acceptance claim from this single-run design. Peak assembly RSS must remain within the 12 GiB allocation; an OOM invalidates execution.

**Stop conditions**

Stop immediately on:

- Binary, input, lockfile, source identity, or QUAST identity mismatch.
- Wrong node, partition, CPU allocation, or memory allocation.
- Existing run outputs or an existing execution lock.
- Input decompression failure.
- Nonzero assembly or QUAST exit.
- Missing/invalid FASTA, stage telemetry, timing receipt, report, or junction record.
- Unexpected subsampling.
- Filesystem errors, exhausted storage, timeout, cancellation, or OOM.

Preserve all partial artifacts. Do not overwrite, resume, silently retry, increase resources, rebuild binaries, or change parameters within this run.

A quality regression is evaluated after all four otherwise valid cells complete, allowing both dataset outcomes to be retained. Failure leaves `assembly_metrics.tsv.partial` and `acceptance.tsv`; these are diagnostic evidence, not a passing benchmark.

**Optional 500 bp QUAST sensitivity analysis**

Disabled for the primary execution. If performed, evaluate all four existing assemblies with:

```bash
quast.py -R "$ref" "$contigs" -o "$quast500_dir" \
  --threads 8 --min-contig 500 --silent
```

Use separate directories under `quast-min500/<dataset>/<version>` and a separate `assembly_metrics_min500.tsv`. Do not rerun assemblies, overwrite primary reports, replace the 200 bp results, or select whichever threshold looks better. The assembler already defaults to a 500 bp output cutoff, so identical results may be expected; verify rather than assume.

**Final verification and evidence closure**

After the Slurm job terminates, obtain accounting from an authorized cluster session:

```bash
sacct -j "$JOB_ID" --parsable2 \
  --format=JobID,JobName,Partition,NodeList,AllocCPUS,ReqMem,State,ExitCode,Elapsed,MaxRSS \
  > "$E/slurm_accounting.tsv"
```

Require the allocation and measured steps to have completed successfully, with exit status `0:0` and no failed, cancelled, timed-out, or OOM steps. If accounting is unavailable or incomplete, final verification remains pending. The script’s `PASS_SELECTIVE` is provisional until this check completes.

Verify:

- Exactly four assembly cells and four primary QUAST reports exist.
- All primary metric rows have the intended dataset/reference/version pairing.
- A/B ingest counts and subsampling settings agree within each dataset.
- Before/after checksums pass.
- Safety and contiguity gates are reported separately.
- Raw evidence and final accounting are retained.
- Artifact checksums verify successfully.

After the job ends, checksum the closed Slurm logs, command trace, status file, accounting receipt, plan, and script in a final evidence manifest. These live logs are intentionally excluded from the script’s earlier manifest.

Report the outcome as `PASS_SELECTIVE`, `REJECTED`, or `BLOCKED/INCOMPLETE`, with explicit reasons. A selective pass does **not** establish global C0-02R acceptance or close the remaining F19 queue-capacity test-sensitivity finding.

Reproduction scope remains **REPRODUCIBLE_ON_ORIGINAL_HPC**; self-contained external reproduction is not established.