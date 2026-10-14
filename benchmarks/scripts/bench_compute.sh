#!/usr/bin/env bash
set -uo pipefail

# Single-run Compute Harness for fastBact Benchmarking Suite (FB-BENCH-001)
# Executed on HPC Compute Node under Slurm. Captures /usr/bin/time -v, hardware context,
# slurm resource accounting, QUAST runs, and produces run.json and metrics.json.

CAMPAIGN=""
DATASET_ID=""
TRACK="E2E"
TOOL=""
REP="1"
THREADS="8"
SEED="101"
OUTDIR=""
R1=""
R2=""
REF="NA"
TRUTH_EVENTS="NA"
IS_NEGATIVE="0"
CUSTOM_CMD=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --campaign) CAMPAIGN="$2"; shift 2 ;;
    --dataset) DATASET_ID="$2"; shift 2 ;;
    --track) TRACK="$2"; shift 2 ;;
    --tool) TOOL="$2"; shift 2 ;;
    --rep) REP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --truth-events) TRUTH_EVENTS="$2"; shift 2 ;;
    --is-negative) IS_NEGATIVE="$2"; shift 2 ;;
    --cmd) CUSTOM_CMD="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [ -z "$OUTDIR" ] || [ -z "$TOOL" ] || [ -z "$DATASET_ID" ]; then
  echo "Usage: bench_compute.sh --campaign <camp> --dataset <id> --tool <tool> --outdir <dir> ..."
  exit 1
fi

mkdir -p "$OUTDIR"
TOOL_DIR="$OUTDIR/tool_out"
mkdir -p "$TOOL_DIR"

# 1. Environment metadata
HOST_NAME=$(hostname)
CPU_MODEL=$(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2 | xargs)
MEM_TOTAL_KB=$(grep -m1 "MemTotal" /proc/meminfo | awk '{print $2}')
KERNEL_VER=$(uname -r)
GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null || echo "UNKNOWN")
GIT_DIRTY=false
if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
  GIT_DIRTY=true
fi

cat <<EOF > "$OUTDIR/environment.json"
{
  "hostname": "${HOST_NAME}",
  "cpu_model": "${CPU_MODEL}",
  "mem_total_kib": ${MEM_TOTAL_KB},
  "kernel": "${KERNEL_VER}",
  "slurm_job_id": "${SLURM_JOB_ID:-LOCAL}",
  "slurm_nodelist": "${SLURM_NODELIST:-LOCAL}",
  "git_commit": "${GIT_COMMIT}",
  "git_dirty": ${GIT_DIRTY}
}
EOF

# 2. Command preparation
if [ -n "$CUSTOM_CMD" ]; then
  EXEC_CMD="$CUSTOM_CMD"
else
  case "$TOOL" in
    fastBact*)
      EXEC_CMD="/hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/target/release/fastBact -1 ${R1} -2 ${R2} -o ${TOOL_DIR} -t ${THREADS}"
      ;;
    spades)
      EXEC_CMD="/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/spades.py --isolate -1 ${R1} -2 ${R2} -o ${TOOL_DIR} -t ${THREADS} -m 32"
      ;;
    shovill)
      EXEC_CMD="/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/shovill --R1 ${R1} --R2 ${R2} --outdir ${TOOL_DIR} --cpus ${THREADS} --ram 12 --force"
      ;;
    skesa)
      EXEC_CMD="bash -c '/hpcfs/fhome/caizhh/.conda/envs/LineaBact/bin/skesa --reads ${R1},${R2} --cores ${THREADS} --memory 32 > ${TOOL_DIR}/contigs.fa'"
      ;;
    *)
      echo "Unknown tool: $TOOL"; exit 1 ;;
  esac
fi

echo "$EXEC_CMD" > "$OUTDIR/command.txt"

# 3. Execution under /usr/bin/time -v
TIME_FILE="$OUTDIR/time.txt"
STDOUT_LOG="$OUTDIR/stdout.log"
STDERR_LOG="$OUTDIR/stderr.log"

START_EPOCH=$(date +%s.%N)
set +e
/usr/bin/time -v sh -c "$EXEC_CMD" > "$STDOUT_LOG" 2> "$STDERR_LOG.raw"
CMD_EXIT=$?
set -e
END_EPOCH=$(date +%s.%N)

# Separate time output from tool stderr
sed -n '/Command being timed:/,$p' "$STDERR_LOG.raw" > "$TIME_FILE"
sed '/Command being timed:/,$d' "$STDERR_LOG.raw" > "$STDERR_LOG"
rm -f "$STDERR_LOG.raw"

# Parse time metrics
WALL_TIME=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "$TIME_FILE" | tr -d ' ')
WALL_S="0"
if [[ "$WALL_TIME" =~ ([0-9]+):([0-9]+)\.([0-9]+) ]]; then
  # m:ss.ss
  MINS=${BASH_REMATCH[1]}
  SECS=${BASH_REMATCH[2]}
  FRAC=${BASH_REMATCH[3]}
  WALL_S=$(echo "$MINS * 60 + $SECS + 0.$FRAC" | bc -l)
elif [[ "$WALL_TIME" =~ ([0-9]+):([0-9]+):([0-9]+) ]]; then
  # h:mm:ss
  HRS=${BASH_REMATCH[1]}
  MINS=${BASH_REMATCH[2]}
  SECS=${BASH_REMATCH[3]}
  WALL_S=$(echo "$HRS * 3600 + $MINS * 60 + $SECS" | bc -l)
elif [[ "$WALL_TIME" =~ ([0-9]+):([0-9]+) ]]; then
  MINS=${BASH_REMATCH[1]}
  SECS=${BASH_REMATCH[2]}
  WALL_S=$(echo "$MINS * 60 + $SECS" | bc -l)
else
  WALL_S=$(echo "$END_EPOCH - $START_EPOCH" | bc -l)
fi

WALL_S=$(echo "$WALL_S" | awk '{if ($0 ~ /^\./) print "0"$0; else if ($0 ~ /^-\./) print "-0"substr($0,2); else print $0}')
USER_S=$(awk -F': ' '/User time \(seconds\)/ {print $2}' "$TIME_FILE" || echo "0")
USER_S=$(echo "$USER_S" | awk '{if ($0 ~ /^\./) print "0"$0; else print $0}')
SYS_S=$(awk -F': ' '/System time \(seconds\)/ {print $2}' "$TIME_FILE" || echo "0")
SYS_S=$(echo "$SYS_S" | awk '{if ($0 ~ /^\./) print "0"$0; else print $0}')
MAX_RSS_KB=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$TIME_FILE" || echo "0")

# Slurm accounting
cat <<EOF > "$OUTDIR/slurm.txt"
SLURM_JOB_ID=${SLURM_JOB_ID:-LOCAL}
SLURM_NODELIST=${SLURM_NODELIST:-LOCAL}
SLURM_CPUS_ON_NODE=${SLURM_CPUS_ON_NODE:-${THREADS}}
EXIT_CODE=${CMD_EXIT}
EOF

# 4. Assess Status
STATUS="PASS"
if [ "$IS_NEGATIVE" = "1" ]; then
  # For negative datasets (N-*), non-zero exit or graceful error log is EXPECTED
  if [ $CMD_EXIT -ne 0 ] || grep -qi -E "error|invalid|corrupt|truncated|mismatch|failed" "$STDERR_LOG" "$STDOUT_LOG"; then
    STATUS="PASS"
  else
    STATUS="FAIL"
  fi
else
  if [ $CMD_EXIT -ne 0 ]; then
    STATUS="FAIL"
  fi
fi

# 5. Locate contigs.fa
ASSEMBLY_FA="$TOOL_DIR/contigs.fa"
if [ ! -f "$ASSEMBLY_FA" ] && [ -f "$TOOL_DIR/contigs.fasta" ]; then
  cp "$TOOL_DIR/contigs.fasta" "$ASSEMBLY_FA"
fi

ASSEMBLY_SHA256="null"
if [ -f "$ASSEMBLY_FA" ]; then
  ASSEMBLY_SHA256="\"$(sha256sum "$ASSEMBLY_FA" | awk '{print $1}')\""
fi

# 6. Evaluation (QUAST & Truth analysis) if status is PASS and reference is available
QUAST_REPORT_FILE="null"
if [ "$STATUS" = "PASS" ] && [ -f "$ASSEMBLY_FA" ] && [ "$REF" != "NA" ] && [ -f "$REF" ]; then
  # QUAST Pass 1: standard 500bp
  rm -rf "$OUTDIR/quast_500"
  /hpcfs/fhome/caizhh/.local/bin/quast.py -R "$REF" "$ASSEMBLY_FA" -o "$OUTDIR/quast_500" --min-contig 500 -t "$THREADS" --silent || true

  # QUAST Pass 2: small-contig / plasmid 200bp
  rm -rf "$OUTDIR/quast_200"
  /hpcfs/fhome/caizhh/.local/bin/quast.py -R "$REF" "$ASSEMBLY_FA" -o "$OUTDIR/quast_200" --min-contig 200 -t "$THREADS" --silent || true

  if [ -f "$OUTDIR/quast_500/report.tsv" ]; then
    QUAST_REPORT_FILE="\"quast_500/report.tsv\""
  fi

  # Run metrics collector & truth evaluator
  python3 /hpcfs/fhome/caizhh/Desktop/03_Tool_Development/01_fastBact/benchmarks/scripts/collect_metrics.py \
    --quast-500 "$OUTDIR/quast_500" \
    --quast-200 "$OUTDIR/quast_200" \
    --assembly "$ASSEMBLY_FA" \
    --reference "$REF" \
    --truth-events "$TRUTH_EVENTS" \
    --out-json "$OUTDIR/metrics.json" || true
else
  # Negative run or failed run: write null metrics
  cat <<EOF > "$OUTDIR/metrics.json"
{
  "schema_version": "1.0",
  "extensive_misassemblies": null,
  "local_misassemblies": null,
  "genome_fraction_pct": null,
  "nga50": null,
  "ng50": null,
  "n50": null,
  "largest_contig": null,
  "total_length": null,
  "contig_count_500bp": null,
  "contig_count_200bp": null,
  "duplication_ratio": null,
  "unaligned_length": null,
  "mismatches_per_100kbp": null,
  "indels_per_100kbp": null,
  "junction_precision": null,
  "junction_recall": null,
  "plasmid_recovery_pct": null,
  "plasmid_circularity_correct": null,
  "plasmid_chimera_count": null
}
EOF
fi

# 7. Write run.json
cat <<EOF > "$OUTDIR/run.json"
{
  "campaign_id": "${CAMPAIGN}",
  "dataset_id": "${DATASET_ID}",
  "track": "${TRACK}",
  "tool": "${TOOL}",
  "tool_version": "${TOOL}",
  "git_commit": "${GIT_COMMIT}",
  "git_dirty": ${GIT_DIRTY},
  "threads": ${THREADS},
  "memory_limit_gib": 32,
  "node_class": "${HOST_NAME} (${CPU_MODEL})",
  "slurm_job_id": "${SLURM_JOB_ID:-LOCAL}",
  "repetition": ${REP},
  "seed": ${SEED},
  "command": $(jq -Rs . < "$OUTDIR/command.txt"),
  "exit_code": ${CMD_EXIT},
  "status": "${STATUS}",
  "wall_s": ${WALL_S:-null},
  "cpu_user_s": ${USER_S:-null},
  "cpu_system_s": ${SYS_S:-null},
  "peak_rss_kib": ${MAX_RSS_KB:-null},
  "scratch_peak_bytes": 0,
  "reference_confidence": "$([ "$REF" != "NA" ] && echo "verified" || echo "unavailable")",
  "assembly_sha256": ${ASSEMBLY_SHA256},
  "quast_metrics_file": ${QUAST_REPORT_FILE},
  "evaluation_version": "QUAST 5.2.0"
}
EOF

echo "[bench_compute] Completed ${TOOL} on ${DATASET_ID} (status=${STATUS}, wall=${WALL_S}s, rss=${MAX_RSS_KB} KiB)"
exit 0
