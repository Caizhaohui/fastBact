#!/usr/bin/env bash
set -uo pipefail

echo "=== Running SPAdes baseline on Real Isolate Cohort R1-001..003 ==="
benchmarks/scripts/run_suite.py --datasets R1-001 R1-002 R1-003 --tools spades
echo "=== SPAdes Real Cohort Complete ==="
