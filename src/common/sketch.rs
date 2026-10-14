//! Streaming cardinality and frequency sketches: HyperLogLog & Count-Min Sketch.
//!
//! - `HyperLogLog`: 14-bit precision (16,384 registers) for unique k-mer cardinality estimation.
//! - `CountMinSketch`: 4-row lock-free sketch with 4-bit saturating counters (packed in `AtomicU8`)
//!   supporting multi-threaded frequency counting up to 15.

use std::sync::atomic::{AtomicU8, Ordering};

// ----------------------------------------------------------------------------
// SplitMix64 Hash Generator
// ----------------------------------------------------------------------------

#[inline(always)]
pub fn splitmix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

// ----------------------------------------------------------------------------
// HyperLogLog (14-bit precision, 16,384 registers)
// ----------------------------------------------------------------------------

pub const HLL_PRECISION: usize = 14;
pub const HLL_REGISTERS: usize = 1 << HLL_PRECISION; // 16,384

/// HyperLogLog cardinality estimator using 14-bit register index and 50-bit zero counting.
#[derive(Clone, Debug)]
pub struct HyperLogLog {
    pub registers: Vec<u8>,
}

impl Default for HyperLogLog {
    fn default() -> Self {
        Self::new()
    }
}

impl HyperLogLog {
    /// Creates a new HyperLogLog sketch with 16,384 8-bit registers initialized to zero.
    pub fn new() -> Self {
        Self {
            registers: vec![0u8; HLL_REGISTERS],
        }
    }

    /// Adds a 64-bit hash item into the sketch.
    #[inline(always)]
    pub fn add_hash(&mut self, h: u64) {
        // High 14 bits determine the register index
        let idx = (h >> (64 - HLL_PRECISION)) as usize;
        // Remaining 50 bits for leading zero run count
        let w = h & ((1u64 << (64 - HLL_PRECISION)) - 1);
        let rho = if w == 0 {
            (64 - HLL_PRECISION + 1) as u8
        } else {
            (w.leading_zeros() as usize - HLL_PRECISION + 1) as u8
        };

        if rho > self.registers[idx] {
            self.registers[idx] = rho;
        }
    }

    /// Merges another HyperLogLog sketch into `self`.
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.registers.len(), other.registers.len());
        for i in 0..HLL_REGISTERS {
            if other.registers[i] > self.registers[i] {
                self.registers[i] = other.registers[i];
            }
        }
    }

    /// Estimates the number of unique elements inserted into the sketch.
    pub fn estimate(&self) -> f64 {
        let m = HLL_REGISTERS as f64;
        // Bias correction constant alpha_m for m = 16384
        let alpha_m = 0.7213 / (1.0 + 1.079 / m);

        let mut sum = 0.0f64;
        let mut zero_registers = 0usize;

        for &val in &self.registers {
            sum += 2.0f64.powi(-(val as i32));
            if val == 0 {
                zero_registers += 1;
            }
        }

        let raw_estimate = alpha_m * m * m / sum;

        // Small range correction (LinearCounting)
        if raw_estimate <= 2.5 * m {
            if zero_registers > 0 {
                m * (m / zero_registers as f64).ln()
            } else {
                raw_estimate
            }
        } else if raw_estimate <= (1.0 / 30.0) * (u64::MAX as f64) {
            raw_estimate
        } else {
            // Large range correction
            let two_to_64 = 18446744073709551616.0f64;
            -two_to_64 * (1.0 - raw_estimate / two_to_64).ln()
        }
    }
}

// ----------------------------------------------------------------------------
// CountMinSketch (4 rows, 4-bit saturating counters, lock-free)
// ----------------------------------------------------------------------------

const CMS_ROWS: usize = 4;
const CMS_SEEDS: [u64; 4] = [
    0x243f6a8885a308d3,
    0x13198a2e03707344,
    0xa4093822299f31d0,
    0x082efa98ec4e6c89,
];

/// Lock-free Count-Min Sketch with 4 rows of 4-bit saturating counters (0..15).
/// Two counters are packed into each `AtomicU8`.
pub struct CountMinSketch {
    pub width: usize,        // Number of counters per row (must be power of 2)
    pub width_mask: usize,   // width - 1
    pub table: Vec<AtomicU8>,// (CMS_ROWS * width) / 2 bytes
}

impl CountMinSketch {
    /// Creates a CountMinSketch with specified width per row (rounded up to power of 2).
    /// Default standard width for bacterial genomes is 4,194,304 (8 MB RAM total).
    pub fn new(mut width: usize) -> Self {
        if !width.is_power_of_two() {
            width = width.next_power_of_two();
        }
        let total_counters = CMS_ROWS * width;
        let total_bytes = total_counters / 2;

        let mut table = Vec::with_capacity(total_bytes);
        for _ in 0..total_bytes {
            table.push(AtomicU8::new(0));
        }

        Self {
            width,
            width_mask: width - 1,
            table,
        }
    }

    /// Standard default 8MB sketch suitable for 1M - 100M k-mers.
    pub fn default_8mb() -> Self {
        Self::new(1 << 22) // 4M counters per row -> 8MB total
    }

    /// Dynamically sizes CountMinSketch according to estimated genome size in bp.
    /// Allocates between 4 MB (width 1<<21) for small/test genomes and 256 MB (width 1<<27) for full genomes.
    pub fn for_genome_size(est_genome_size: u64) -> Self {
        let target = (est_genome_size * 16).next_power_of_two() as usize;
        let width = target.clamp(1 << 21, 1 << 27);
        Self::new(width)
    }

    /// Increments the counter for a given 64-bit key across all 4 rows.
    /// Counters saturate at 15. Lock-free and safe across multiple rayon threads.
    #[inline(always)]
    pub fn add(&self, key: u64) {
        for r in 0..CMS_ROWS {
            let h = splitmix64(key ^ CMS_SEEDS[r]);
            let col = (h as usize) & self.width_mask;
            let counter_idx = r * self.width + col;
            let byte_idx = counter_idx >> 1;
            let is_high = (counter_idx & 1) == 1;

            let atomic_byte = &self.table[byte_idx];
            let mut cur = atomic_byte.load(Ordering::Relaxed);
            loop {
                let cnt = if is_high { (cur >> 4) & 0x0F } else { cur & 0x0F };
                if cnt >= 15 {
                    break;
                }
                let next_cnt = cnt + 1;
                let next_val = if is_high {
                    (cur & 0x0F) | (next_cnt << 4)
                } else {
                    (cur & 0xF0) | next_cnt
                };
                match atomic_byte.compare_exchange_weak(
                    cur,
                    next_val,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual) => cur = actual,
                }
            }
        }
    }

    /// Estimates the frequency count of `key`. Returns a value in 0..15.
    #[inline(always)]
    pub fn estimate(&self, key: u64) -> u8 {
        let mut min_val = 15u8;
        for r in 0..CMS_ROWS {
            let h = splitmix64(key ^ CMS_SEEDS[r]);
            let col = (h as usize) & self.width_mask;
            let counter_idx = r * self.width + col;
            let byte_idx = counter_idx >> 1;
            let is_high = (counter_idx & 1) == 1;

            let cur = self.table[byte_idx].load(Ordering::Relaxed);
            let cnt = if is_high { (cur >> 4) & 0x0F } else { cur & 0x0F };
            if cnt < min_val {
                min_val = cnt;
            }
            if min_val == 0 {
                return 0;
            }
        }
        min_val
    }

    /// Increments and returns the updated estimated frequency.
    #[inline(always)]
    pub fn add_and_estimate(&self, key: u64) -> u8 {
        self.add(key);
        self.estimate(key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn test_hll_cardinality_verification_gate() {
        // Verification Gate item 3:
        // HyperLogLog estimates unique 31-mers within +/- 3% error on 100,000 synthetic k-mers.
        let mut hll = HyperLogLog::new();
        let target_unique: usize = 100_000;

        for i in 0..target_unique {
            // Generate synthetic 31-mer 64-bit hash
            let h = splitmix64(i as u64);
            hll.add_hash(h);
        }

        let estimated = hll.estimate();
        let error = (estimated - target_unique as f64).abs() / (target_unique as f64);

        assert!(
            error < 0.03,
            "HLL error too high: estimated={}, target={}, error={:.4} (must be < 0.03)",
            estimated,
            target_unique,
            error
        );
    }

    #[test]
    fn test_hll_merge() {
        let mut hll1 = HyperLogLog::new();
        let mut hll2 = HyperLogLog::new();

        for i in 0..50_000 {
            hll1.add_hash(splitmix64(i as u64));
        }
        for i in 50_000..100_000 {
            hll2.add_hash(splitmix64(i as u64));
        }

        hll1.merge(&hll2);
        let estimated = hll1.estimate();
        let error = (estimated - 100_000.0).abs() / 100_000.0;
        assert!(error < 0.03, "Merged HLL error too high: {}", error);
    }

    #[test]
    fn test_count_min_sketch_basic() {
        let cms = CountMinSketch::new(1024);
        assert_eq!(cms.estimate(42), 0);

        for _ in 0..5 {
            cms.add(42);
        }
        assert_eq!(cms.estimate(42), 5);

        for _ in 0..20 {
            cms.add(42);
        }
        // Counter saturates at 15
        assert_eq!(cms.estimate(42), 15);
    }

    #[test]
    fn test_count_min_sketch_multithreaded() {
        let cms = Arc::new(CountMinSketch::new(4096));
        let num_threads = 4;
        let iters_per_thread = 3;
        let mut handles = Vec::new();

        for _ in 0..num_threads {
            let cms_clone = Arc::clone(&cms);
            handles.push(std::thread::spawn(move || {
                for _ in 0..iters_per_thread {
                    cms_clone.add(12345);
                }
            }));
        }

        for h in handles {
            h.join().unwrap();
        }

        assert_eq!(cms.estimate(12345), (num_threads * iters_per_thread) as u8);
    }
}
