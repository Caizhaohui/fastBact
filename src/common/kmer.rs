//! Compile-time polymorphic K-mer representations, rolling hashes, and minimizers.
//!
//! Provides the `KmerInt` trait with implementations for:
//! - `u64`: k <= 31
//! - `u128`: 33 <= k <= 63
//! - `[u64; 4]`: 65 <= k <= 127
//!
//! Features O(1) rolling canonical k-mer updates and `ntHash` 64-bit rolling hash.

use crate::common::bitseq::{
    rev_comp_word, ASCII_TO_2BIT, COMPLEMENT_2BIT, TWO_BIT_TO_ASCII,
};
pub use crate::common::sketch::splitmix64;
use std::fmt::Debug;
use std::hash::Hash;

/// Trait for compact integer representations of DNA k-mers.
pub trait KmerInt:
    Copy + Clone + Eq + Ord + Hash + Send + Sync + 'static + Debug
{
    /// Maximum k-mer size supported by this integer type.
    const MAX_K: usize;

    /// Zero representation.
    fn zero() -> Self;

    /// Encodes an ASCII DNA string into a k-mer.
    fn from_ascii(ascii: &[u8]) -> Self;

    /// Decodes the k-mer into an ASCII DNA string of length `k`.
    fn to_ascii(&self, k: usize) -> Vec<u8>;

    /// Appends a 2-bit base (0..3) at the 3' end, shifting existing bases left and masking to `k`.
    fn push_base(&self, base: u8, k: usize) -> Self;

    /// Returns the reverse complement of this k-mer of length `k`.
    fn rev_comp(&self, k: usize) -> Self;

    /// Returns the canonical k-mer (`min(kmer, rc)`) and a boolean flag indicating
    /// whether the reverse complement was chosen (`true` if rc < fwd).
    #[inline(always)]
    fn canonical(&self, k: usize) -> (Self, bool) {
        let rc = self.rev_comp(k);
        if rc < *self {
            (rc, true)
        } else {
            (*self, false)
        }
    }

    /// Performs an O(1) rolling update for both forward and reverse complement simultaneously.
    ///
    /// Given the previous `(fwd, rev)` pair and the incoming 3' base (0..3),
    /// produces the next `(fwd', rev')` in O(1) without iterating over all bases.
    fn roll_next(fwd: Self, rev: Self, base: u8, k: usize) -> (Self, Self);

    /// Converts or hashes the k-mer to a 64-bit unsigned integer.
    fn to_u64(&self) -> u64;
}

// ----------------------------------------------------------------------------
// u64 Implementation (k <= 31)
// ----------------------------------------------------------------------------

impl KmerInt for u64 {
    const MAX_K: usize = 31;

    #[inline(always)]
    fn zero() -> Self {
        0
    }

    fn from_ascii(ascii: &[u8]) -> Self {
        let k = ascii.len();
        debug_assert!(k <= Self::MAX_K);
        let mut val = 0u64;
        for &b in ascii {
            val = (val << 2) | (ASCII_TO_2BIT[b as usize] as u64);
        }
        val
    }

    fn to_ascii(&self, k: usize) -> Vec<u8> {
        let mut buf = Vec::with_capacity(k);
        for i in (0..k).rev() {
            let code = ((*self >> (2 * i)) & 0b11) as usize;
            buf.push(TWO_BIT_TO_ASCII[code]);
        }
        buf
    }

    #[inline(always)]
    fn push_base(&self, base: u8, k: usize) -> Self {
        let mask = if k >= 32 { !0u64 } else { (1u64 << (2 * k)) - 1 };
        ((*self << 2) | ((base & 0b11) as u64)) & mask
    }

    #[inline(always)]
    fn rev_comp(&self, k: usize) -> Self {
        rev_comp_word(*self, k)
    }

    #[inline(always)]
    fn roll_next(fwd: Self, rev: Self, base: u8, k: usize) -> (Self, Self) {
        let b = (base & 0b11) as u64;
        let comp_b = (COMPLEMENT_2BIT[b as usize] as u64) << (2 * (k - 1));
        let mask = if k >= 32 { !0u64 } else { (1u64 << (2 * k)) - 1 };

        let next_fwd = ((fwd << 2) | b) & mask;
        let next_rev = (rev >> 2) | comp_b;
        (next_fwd, next_rev)
    }

    #[inline(always)]
    fn to_u64(&self) -> u64 {
        *self
    }
}

// ----------------------------------------------------------------------------
// u128 Implementation (33 <= k <= 63)
// ----------------------------------------------------------------------------

impl KmerInt for u128 {
    const MAX_K: usize = 63;

    #[inline(always)]
    fn zero() -> Self {
        0
    }

    fn from_ascii(ascii: &[u8]) -> Self {
        let k = ascii.len();
        debug_assert!(k <= Self::MAX_K);
        let mut val = 0u128;
        for &b in ascii {
            val = (val << 2) | (ASCII_TO_2BIT[b as usize] as u128);
        }
        val
    }

    fn to_ascii(&self, k: usize) -> Vec<u8> {
        let mut buf = Vec::with_capacity(k);
        for i in (0..k).rev() {
            let code = ((*self >> (2 * i)) & 0b11) as usize;
            buf.push(TWO_BIT_TO_ASCII[code]);
        }
        buf
    }

    #[inline(always)]
    fn push_base(&self, base: u8, k: usize) -> Self {
        let mask = if k >= 64 { !0u128 } else { (1u128 << (2 * k)) - 1 };
        ((*self << 2) | ((base & 0b11) as u128)) & mask
    }

    #[inline(always)]
    fn rev_comp(&self, k: usize) -> Self {
        if k == 0 {
            return 0;
        }
        let hi = (*self >> 64) as u64;
        let lo = *self as u64;

        let rc_hi = rev_comp_word(hi, 32);
        let rc_lo = rev_comp_word(lo, 32);

        // Swap 64-bit halves
        let full_rc = ((rc_lo as u128) << 64) | (rc_hi as u128);

        if k < 64 {
            let shift = (64 - k) * 2;
            let mask = (1u128 << (2 * k)) - 1;
            (full_rc >> shift) & mask
        } else {
            full_rc
        }
    }

    #[inline(always)]
    fn roll_next(fwd: Self, rev: Self, base: u8, k: usize) -> (Self, Self) {
        let b = (base & 0b11) as u128;
        let comp_b = (COMPLEMENT_2BIT[b as usize] as u128) << (2 * (k - 1));
        let mask = if k >= 64 { !0u128 } else { (1u128 << (2 * k)) - 1 };

        let next_fwd = ((fwd << 2) | b) & mask;
        let next_rev = (rev >> 2) | comp_b;
        (next_fwd, next_rev)
    }

    #[inline(always)]
    fn to_u64(&self) -> u64 {
        (*self as u64) ^ ((*self >> 64) as u64)
    }
}

// ----------------------------------------------------------------------------
// [u64; 4] Implementation (65 <= k <= 127)
// Stored big-endian: words[0] contains MSB (bases 96..127), words[3] contains LSB (bases 0..31).
// ----------------------------------------------------------------------------

impl KmerInt for [u64; 4] {
    const MAX_K: usize = 127;

    #[inline(always)]
    fn zero() -> Self {
        [0, 0, 0, 0]
    }

    fn from_ascii(ascii: &[u8]) -> Self {
        let k = ascii.len();
        debug_assert!(k <= Self::MAX_K);
        let mut words = [0u64; 4];
        for &b in ascii {
            let code = ASCII_TO_2BIT[b as usize] as u64;
            words[0] = (words[0] << 2) | (words[1] >> 62);
            words[1] = (words[1] << 2) | (words[2] >> 62);
            words[2] = (words[2] << 2) | (words[3] >> 62);
            words[3] = (words[3] << 2) | code;
        }
        words
    }

    fn to_ascii(&self, k: usize) -> Vec<u8> {
        let mut buf = Vec::with_capacity(k);
        for i in (0..k).rev() {
            let word_from_lsb = i / 32;
            let word_idx = 3 - word_from_lsb;
            let shift = (i % 32) * 2;
            let code = ((self[word_idx] >> shift) & 0b11) as usize;
            buf.push(TWO_BIT_TO_ASCII[code]);
        }
        buf
    }

    fn push_base(&self, base: u8, k: usize) -> Self {
        let code = (base & 0b11) as u64;
        let mut words = [
            (self[0] << 2) | (self[1] >> 62),
            (self[1] << 2) | (self[2] >> 62),
            (self[2] << 2) | (self[3] >> 62),
            (self[3] << 2) | code,
        ];

        // Mask to 2*k bits
        let total_bits = 2 * k;
        for j in 0..4 {
            let word_min_bit = (3 - j) * 64;
            if total_bits <= word_min_bit {
                words[j] = 0;
            } else if total_bits < word_min_bit + 64 {
                let bits_in_word = total_bits - word_min_bit;
                let mask = (1u64 << bits_in_word) - 1;
                words[j] &= mask;
            }
        }
        words
    }

    fn rev_comp(&self, k: usize) -> Self {
        if k == 0 {
            return [0, 0, 0, 0];
        }
        // Reverse-complement each 32-base word and swap order
        let full_rc = [
            rev_comp_word(self[3], 32),
            rev_comp_word(self[2], 32),
            rev_comp_word(self[1], 32),
            rev_comp_word(self[0], 32),
        ];

        if k == 128 {
            return full_rc;
        }

        // Shift right by 2 * (128 - k)
        let s = (128 - k) * 2;
        let w_shift = s / 64;
        let b_shift = s % 64;

        let mut shifted = [0u64; 4];
        for i in 0..4 {
            if i + w_shift < 4 {
                shifted[i + w_shift] |= full_rc[i] >> b_shift;
                if b_shift > 0 && i + w_shift + 1 < 4 {
                    shifted[i + w_shift + 1] |= full_rc[i] << (64 - b_shift);
                }
            }
        }

        // Mask top unused bits
        let total_bits = 2 * k;
        for j in 0..4 {
            let word_min_bit = (3 - j) * 64;
            if total_bits <= word_min_bit {
                shifted[j] = 0;
            } else if total_bits < word_min_bit + 64 {
                let bits_in_word = total_bits - word_min_bit;
                let mask = (1u64 << bits_in_word) - 1;
                shifted[j] &= mask;
            }
        }

        shifted
    }

    fn roll_next(fwd: Self, rev: Self, base: u8, k: usize) -> (Self, Self) {
        let next_fwd = fwd.push_base(base, k);

        let b = (base & 0b11) as u64;
        let comp_b = COMPLEMENT_2BIT[b as usize] as u64;

        // Shift rev right by 2
        let mut next_rev = [
            rev[0] >> 2,
            (rev[1] >> 2) | (rev[0] << 62),
            (rev[2] >> 2) | (rev[1] << 62),
            (rev[3] >> 2) | (rev[2] << 62),
        ];

        // Insert comp_b at bit position 2 * (k - 1)
        let bit_pos = 2 * (k - 1);
        let word_from_lsb = bit_pos / 64;
        let word_idx = 3 - word_from_lsb;
        let shift = bit_pos % 64;
        next_rev[word_idx] |= comp_b << shift;

        (next_fwd, next_rev)
    }

    #[inline(always)]
    fn to_u64(&self) -> u64 {
        self[0] ^ self[1] ^ self[2] ^ self[3]
    }
}

// ----------------------------------------------------------------------------
// RollingKmer Struct
// ----------------------------------------------------------------------------

/// Helper struct that maintains both forward and reverse-complement orientations
/// in O(1) time as new bases are streamed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RollingKmer<K: KmerInt> {
    pub k: usize,
    pub fwd: K,
    pub rev: K,
}

impl<K: KmerInt> RollingKmer<K> {
    /// Initializes a rolling k-mer from an initial ASCII sequence of length `k`.
    pub fn from_ascii(ascii: &[u8]) -> Self {
        let k = ascii.len();
        let fwd = K::from_ascii(ascii);
        let rev = fwd.rev_comp(k);
        Self { k, fwd, rev }
    }

    /// Advances the window by one base, updating `fwd` and `rev` in O(1).
    /// Returns `(canonical_kmer, is_reverse_complement)`.
    #[inline(always)]
    pub fn roll_next(&mut self, base: u8) -> (K, bool) {
        let (fwd, rev) = K::roll_next(self.fwd, self.rev, base, self.k);
        self.fwd = fwd;
        self.rev = rev;
        self.canonical()
    }

    /// Returns the canonical k-mer and orientation flag.
    #[inline(always)]
    pub fn canonical(&self) -> (K, bool) {
        if self.rev < self.fwd {
            (self.rev, true)
        } else {
            (self.fwd, false)
        }
    }
}

// ----------------------------------------------------------------------------
// ntHash: 64-bit Rolling Hash for DNA Sequences
// ----------------------------------------------------------------------------

/// 64-bit random values for DNA bases A, C, G, T (from ntHash specification).
pub const NT_HASH_A: u64 = 0x3c8bfbb395c60474;
pub const NT_HASH_C: u64 = 0x3193c18562a02b4c;
pub const NT_HASH_G: u64 = 0x20323ee046061444;
pub const NT_HASH_T: u64 = 0x295c504fac61002e;

pub static NT_HASH_FWD: [u64; 4] = [NT_HASH_A, NT_HASH_C, NT_HASH_G, NT_HASH_T];
pub static NT_HASH_REV: [u64; 4] = [NT_HASH_T, NT_HASH_G, NT_HASH_C, NT_HASH_A];

/// Computes the initial 64-bit forward ntHash for a DNA slice of length `k`.
pub fn nthash_fwd(seq: &[u8]) -> u64 {
    let k = seq.len();
    let mut h = 0u64;
    for (i, &b) in seq.iter().enumerate() {
        let code = ASCII_TO_2BIT[b as usize] as usize;
        let val = NT_HASH_FWD[code];
        let shift = (k - 1 - i) % 64;
        h ^= val.rotate_left(shift as u32);
    }
    h
}

/// Computes the initial 64-bit reverse-complement ntHash for a DNA slice of length `k`.
pub fn nthash_rev(seq: &[u8]) -> u64 {
    let mut h = 0u64;
    for (i, &b) in seq.iter().enumerate() {
        let code = ASCII_TO_2BIT[b as usize] as usize;
        let val = NT_HASH_REV[code];
        let shift = i % 64;
        h ^= val.rotate_left(shift as u32);
    }
    h
}

/// Computes the canonical 64-bit ntHash for a sequence: `min(fwd, rev)`.
pub fn nthash_canonical(seq: &[u8]) -> u64 {
    let f = nthash_fwd(seq);
    let r = nthash_rev(seq);
    f.min(r)
}

/// Stateful rolling ntHash iterator.
#[derive(Clone, Copy, Debug)]
pub struct RollingNtHash {
    pub k: usize,
    pub h_fwd: u64,
    pub h_rev: u64,
}

impl RollingNtHash {
    /// Initializes rolling hash from initial k-mer ASCII slice.
    pub fn new(seq: &[u8]) -> Self {
        Self {
            k: seq.len(),
            h_fwd: nthash_fwd(seq),
            h_rev: nthash_rev(seq),
        }
    }

    /// Rolls the hash window by one base in O(1).
    /// `out_base`: base dropping out of the 5' end.
    /// `in_base`: base entering at the 3' end.
    #[inline(always)]
    pub fn roll(&mut self, out_base: u8, in_base: u8) -> u64 {
        let out_code = ASCII_TO_2BIT[out_base as usize] as usize;
        let in_code = ASCII_TO_2BIT[in_base as usize] as usize;

        // Forward roll: rol(h_fwd, 1) ^ rol(H(out), k % 64) ^ H(in)
        let k_mod = (self.k % 64) as u32;
        self.h_fwd = self.h_fwd.rotate_left(1)
            ^ NT_HASH_FWD[out_code].rotate_left(k_mod)
            ^ NT_HASH_FWD[in_code];

        // Reverse roll: ror(h_rev, 1) ^ ror(H(rc(out)), 1) ^ rol(H(rc(in)), (k - 1) % 64)
        let km1_mod = ((self.k - 1) % 64) as u32;
        self.h_rev = self.h_rev.rotate_right(1)
            ^ NT_HASH_REV[out_code].rotate_right(1)
            ^ NT_HASH_REV[in_code].rotate_left(km1_mod);

        self.canonical()
    }

    /// Returns canonical hash `min(h_fwd, h_rev)`.
    #[inline(always)]
    pub fn canonical(&self) -> u64 {
        self.h_fwd.min(self.h_rev)
    }

    /// Returns canonical hash mixed via SplitMix64 for unbiased uniform distribution.
    #[inline(always)]
    pub fn canonical_mixed(&self) -> u64 {
        splitmix64(self.canonical())
    }
}

// ----------------------------------------------------------------------------
// Sliding Window Minimizers
// ----------------------------------------------------------------------------

/// Extracts canonical minimizer hashes from a DNA sequence with minimizer size `m` (default 15).
/// Returns a vector of `(minimizer_hash, start_pos)`.
pub fn extract_minimizers(seq: &[u8], m: usize, w: usize) -> Vec<(u64, usize)> {
    if seq.len() < m {
        return Vec::new();
    }

    let num_m_mers = seq.len() - m + 1;
    let mut m_hashes = Vec::with_capacity(num_m_mers);

    // Compute all m-mer hashes
    let mut rolling = RollingNtHash::new(&seq[..m]);
    m_hashes.push((rolling.canonical(), 0));

    for i in 1..num_m_mers {
        let h = rolling.roll(seq[i - 1], seq[i + m - 1]);
        m_hashes.push((h, i));
    }

    if m_hashes.len() <= w {
        // Return single minimum across the whole sequence
        if let Some(&min_entry) = m_hashes.iter().min_by_key(|(h, _)| *h) {
            return vec![min_entry];
        }
        return Vec::new();
    }

    let mut minimizers = Vec::new();
    let mut min_val = u64::MAX;
    let mut min_idx = 0;

    for i in 0..m_hashes.len() {
        if i >= w && min_idx <= i - w {
            // Min fell out of window, recompute in [i - w + 1 ..= i]
            min_val = u64::MAX;
            for j in (i + 1 - w)..=i {
                if m_hashes[j].0 < min_val {
                    min_val = m_hashes[j].0;
                    min_idx = j;
                }
            }
            minimizers.push(m_hashes[min_idx]);
        } else {
            if m_hashes[i].0 < min_val {
                min_val = m_hashes[i].0;
                min_idx = i;
                minimizers.push(m_hashes[min_idx]);
            }
        }
    }

    minimizers.dedup_by(|a, b| a.1 == b.1);
    minimizers
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::bitseq::reverse_complement_ascii;

    #[test]
    fn test_u64_kmer_basic() {
        let seq = b"ACGTACGTACGTACGTACGTACGT"; // 24 bp
        let k = 24;
        let kmer = <u64 as KmerInt>::from_ascii(seq);
        let unpacked = kmer.to_ascii(k);
        assert_eq!(unpacked, seq);

        let rc = kmer.rev_comp(k);
        assert_eq!(rc.to_ascii(k), reverse_complement_ascii(seq));
    }

    #[test]
    fn test_u128_kmer_basic() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 56 bp
        let k = 56;
        let kmer = <u128 as KmerInt>::from_ascii(seq);
        let unpacked = kmer.to_ascii(k);
        assert_eq!(unpacked, seq);

        let rc = kmer.rev_comp(k);
        assert_eq!(rc.to_ascii(k), reverse_complement_ascii(seq));
    }

    #[test]
    fn test_u256_kmer_basic() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 96 bp
        let k = 96;
        let kmer = <[u64; 4]>::from_ascii(seq);
        let unpacked = kmer.to_ascii(k);
        assert_eq!(unpacked, seq);

        let rc = kmer.rev_comp(k);
        assert_eq!(rc.to_ascii(k), reverse_complement_ascii(seq));
    }

    #[test]
    fn test_rolling_kmer_verification_gate() {
        // Verification Gate item 2:
        // Rolling canonical k-mer extraction for k=21, 55, 99, 127 yields identical
        // canonical representations for a sequence and its reverse complement.

        fn verify_for_k<K: KmerInt>(seq: &[u8], k: usize) {
            let rc_seq = reverse_complement_ascii(seq);
            let n = seq.len();

            // Extract canonical k-mers along forward sequence
            let mut fwd_rolling = RollingKmer::<K>::from_ascii(&seq[..k]);
            let mut fwd_canonicals = Vec::new();
            fwd_canonicals.push(fwd_rolling.canonical().0);

            for i in k..n {
                let (c, _) = fwd_rolling.roll_next(ASCII_TO_2BIT[seq[i] as usize]);
                fwd_canonicals.push(c);
            }

            // Extract canonical k-mers along reverse complement sequence
            let mut rc_rolling = RollingKmer::<K>::from_ascii(&rc_seq[..k]);
            let mut rc_canonicals = Vec::new();
            rc_canonicals.push(rc_rolling.canonical().0);

            for i in k..n {
                let (c, _) = rc_rolling.roll_next(ASCII_TO_2BIT[rc_seq[i] as usize]);
                rc_canonicals.push(c);
            }

            // Reverse the RC list to match coordinates
            rc_canonicals.reverse();

            assert_eq!(
                fwd_canonicals.len(),
                rc_canonicals.len(),
                "Length mismatch for k={}",
                k
            );

            for i in 0..fwd_canonicals.len() {
                assert_eq!(
                    fwd_canonicals[i], rc_canonicals[i],
                    "Mismatch at window {} for k={}",
                    i, k
                );
            }
        }

        let alphabet = [b'A', b'C', b'G', b'T'];
        let mut rng: u64 = 0xCAFE_BABE_DEAD_BEEF;
        let mut rand_base = || {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            alphabet[(rng % 4) as usize]
        };

        let mut long_seq = Vec::with_capacity(300);
        for _ in 0..300 {
            long_seq.push(rand_base());
        }

        // Test k = 21 with u64
        verify_for_k::<u64>(&long_seq, 21);

        // Test k = 55 with u128
        verify_for_k::<u128>(&long_seq, 55);

        // Test k = 99 with [u64; 4]
        verify_for_k::<[u64; 4]>(&long_seq, 99);

        // Test k = 127 with [u64; 4]
        verify_for_k::<[u64; 4]>(&long_seq, 127);
    }

    #[test]
    fn test_nthash_rolling() {
        let seq = b"ACGTACGTACGTTTGCAACCGGTTACGT";
        let k = 15;
        let mut rolling = RollingNtHash::new(&seq[..k]);

        for i in 1..=(seq.len() - k) {
            let roll_h = rolling.roll(seq[i - 1], seq[i + k - 1]);
            let fresh_h = nthash_canonical(&seq[i..i + k]);
            assert_eq!(roll_h, fresh_h, "ntHash mismatch at offset {}", i);
        }
    }

    #[test]
    fn test_minimizers_extraction() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let minimizers = extract_minimizers(seq, 15, 10);
        assert!(!minimizers.is_empty());
    }
}
