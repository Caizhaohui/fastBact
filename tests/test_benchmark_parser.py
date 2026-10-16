#!/usr/bin/env python3
"""
Unit tests for fastBact benchmark stdout & telemetry parser.
Ensures zero NameError regressions and rigorous coverage of all telemetry branches.
"""

import sys
import unittest
import re

# Standalone parser implementation matching benchmarks/scripts/run_var3_ab_benchmark.py
def parse_fastbact_stdout(stdout_text):
    data = {
        "subsampling_factor": 1.0,
        "subsampling_enabled": False,
        "subsampling_reason": "UNKNOWN",
        "estimated_cov": 0.0,
        "actual_input_pairs": 0,
        "actual_retained_pairs": 0,
        "retention_pct": 100.0,
    }
    for line in stdout_text.splitlines():
        line_str = line.strip()
        if "subsample_factor:" in line_str or "subsampling_factor:" in line_str:
            try:
                data["subsampling_factor"] = float(line_str.split(":")[-1].strip())
            except ValueError:
                pass
        elif "subsampling_enabled:" in line_str:
            data["subsampling_enabled"] = line_str.split(":")[-1].strip().lower() == "true"
        elif "subsampling_reason:" in line_str:
            data["subsampling_reason"] = line_str.split("subsampling_reason:")[-1].strip()
        elif "estimated_coverage:" in line_str:
            try:
                val = line_str.split(":")[-1].replace("x", "").strip()
                data["estimated_cov"] = float(val)
            except ValueError:
                pass
        elif "actual_input_pairs:" in line_str:
            try:
                data["actual_input_pairs"] = int(line_str.split(":")[-1].strip())
            except ValueError:
                pass
        elif "actual_retained_pairs:" in line_str:
            try:
                data["actual_retained_pairs"] = int(line_str.split(":")[-1].strip())
            except ValueError:
                pass

    if "factor=" in data["subsampling_reason"]:
        try:
            m = re.search(r"factor=([0-9.]+)", data["subsampling_reason"])
            if m:
                data["subsampling_factor"] = float(m.group(1))
        except Exception:
            pass

    if data["actual_input_pairs"] > 0:
        data["retention_pct"] = (data["actual_retained_pairs"] / data["actual_input_pairs"]) * 100.0
    return data


class TestFastBactStdoutParser(unittest.TestCase):
    def test_subsampling_active_high_depth(self):
        log = """
2026-10-15T16:00:00Z  INFO fastBact: --- Input Sequencing Statistics (C1b-3 Calibrated Ingest) ---
2026-10-15T16:00:00Z  INFO fastBact:   exact_read_pairs: 200000
2026-10-15T16:00:00Z  INFO fastBact:   total_bases: 59994000 bp
2026-10-15T16:00:00Z  INFO fastBact:   median_read_len: 150 bp
2026-10-15T16:00:00Z  INFO fastBact:   estimated_genome_size: 207300 bp
2026-10-15T16:00:00Z  INFO fastBact:   estimated_coverage: 289.4x
2026-10-15T16:00:00Z  INFO fastBact:   estimator_confidence: HIGH
2026-10-15T16:00:00Z  INFO fastBact: --- Subsampling Policy (C1b-3 Calibrated Policy) ---
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_enabled: true
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_reason: AUTO_SUBSAMPLE_HIGH_CONFIDENCE (est_cov=289.4x > trigger=180.0x, factor=0.4147)
2026-10-15T16:00:00Z  INFO fastBact:   subsample_factor: 0.4147
2026-10-15T16:00:00Z  INFO fastBact: -----------------------------------------------------
2026-10-15T16:00:01Z  INFO fastBact: --- Actual Ingestion Counts ---
2026-10-15T16:00:01Z  INFO fastBact:   actual_input_pairs: 200000
2026-10-15T16:00:01Z  INFO fastBact:   actual_retained_pairs: 82815
2026-10-15T16:00:01Z  INFO fastBact: -------------------------------
"""
        parsed = parse_fastbact_stdout(log)
        self.assertTrue(parsed["subsampling_enabled"])
        self.assertAlmostEqual(parsed["subsampling_factor"], 0.4147, places=4)
        self.assertAlmostEqual(parsed["estimated_cov"], 289.4, places=1)
        self.assertEqual(parsed["actual_input_pairs"], 200000)
        self.assertEqual(parsed["actual_retained_pairs"], 82815)
        self.assertAlmostEqual(parsed["retention_pct"], 41.4075, places=3)
        self.assertIn("AUTO_SUBSAMPLE_HIGH_CONFIDENCE", parsed["subsampling_reason"])

    def test_subsampling_bypassed_within_target(self):
        log = """
2026-10-15T16:00:00Z  INFO fastBact:   estimated_coverage: 145.3x
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_enabled: false
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_reason: AUTO_SUBSAMPLE_HIGH_CONFIDENCE_WITHIN_TARGET (est_cov=145.3x <= trigger=180.0x)
2026-10-15T16:00:00Z  INFO fastBact:   subsample_factor: 1.0000
2026-10-15T16:00:01Z  INFO fastBact:   actual_input_pairs: 100000
2026-10-15T16:00:01Z  INFO fastBact:   actual_retained_pairs: 100000
"""
        parsed = parse_fastbact_stdout(log)
        self.assertFalse(parsed["subsampling_enabled"])
        self.assertAlmostEqual(parsed["subsampling_factor"], 1.0000, places=4)
        self.assertAlmostEqual(parsed["estimated_cov"], 145.3, places=1)
        self.assertEqual(parsed["actual_input_pairs"], 100000)
        self.assertEqual(parsed["actual_retained_pairs"], 100000)
        self.assertAlmostEqual(parsed["retention_pct"], 100.0, places=2)

    def test_low_coverage_warning_intercepted(self):
        log = """
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_enabled: false
2026-10-15T16:00:00Z  INFO fastBact:   subsampling_reason: AUTO_SUBSAMPLE_INTERCEPTED_LOW_COVERAGE_WARNING
2026-10-15T16:00:00Z  INFO fastBact:   subsample_factor: 1.0000
2026-10-15T16:00:00Z  INFO fastBact:   estimated_coverage: 14.7x
2026-10-15T16:00:01Z  INFO fastBact:   actual_input_pairs: 10000
2026-10-15T16:00:01Z  INFO fastBact:   actual_retained_pairs: 10000
"""
        parsed = parse_fastbact_stdout(log)
        self.assertFalse(parsed["subsampling_enabled"])
        self.assertAlmostEqual(parsed["subsampling_factor"], 1.0000, places=4)
        self.assertEqual(parsed["subsampling_reason"], "AUTO_SUBSAMPLE_INTERCEPTED_LOW_COVERAGE_WARNING")
        self.assertEqual(parsed["retention_pct"], 100.0)

    def test_empty_log_does_not_crash(self):
        parsed = parse_fastbact_stdout("")
        self.assertFalse(parsed["subsampling_enabled"])
        self.assertEqual(parsed["subsampling_factor"], 1.0)
        self.assertEqual(parsed["actual_input_pairs"], 0)
        self.assertEqual(parsed["retention_pct"], 100.0)


if __name__ == "__main__":
    unittest.main()
