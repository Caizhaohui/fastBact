//! 2-bit DNA sequence encoding and compact arena storage.
//!
//! Encodes DNA bases into 2 bits:
//! - A = 0b00 (0)
//! - C = 0b01 (1)
//! - G = 0b10 (2)
//! - T = 0b11 (3)
//!
//! Compacts sequences into `Vec<u64>` (32 bases per word) with zero per-read
//! heap allocations.

/// 2-bit DNA encoding: A = 0b00, C = 0b01, G = 0b10, T = 0b11
#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    #[inline(always)]
    pub fn from_ascii(b: u8) -> Self {
        match b {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            _ => Base::A, // Deterministic mapping for N and ambiguous bases
        }
    }

    #[inline(always)]
    pub fn to_ascii(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
        }
    }

    #[inline(always)]
    pub fn complement(self) -> Self {
        match self {
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
        }
    }

    #[inline(always)]
    pub fn from_2bit(code: u8) -> Self {
        match code & 0b11 {
            0 => Base::A,
            1 => Base::C,
            2 => Base::G,
            3 => Base::T,
            _ => unreachable!(),
        }
    }
}

/// Lookup table for mapping any ASCII byte directly to its 2-bit code (0..3).
pub static ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        table[i] = match i as u8 {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0, // Deterministic mapping for non-ACGT bases
        };
        i += 1;
    }
    table
};

/// 2-bit complement lookup table: 0->3, 1->2, 2->1, 3->0.
pub static COMPLEMENT_2BIT: [u8; 4] = [3, 2, 1, 0];

/// 2-bit code to uppercase ASCII base.
pub static TWO_BIT_TO_ASCII: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Reverse-complements an ASCII DNA slice, allocating a new `Vec<u8>`.
pub fn reverse_complement_ascii(ascii: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(ascii.len());
    for &b in ascii.iter().rev() {
        let code = ASCII_TO_2BIT[b as usize];
        let comp_code = COMPLEMENT_2BIT[code as usize];
        rc.push(TWO_BIT_TO_ASCII[comp_code as usize]);
    }
    rc
}

/// Reverses and complements an ASCII DNA slice in-place.
pub fn reverse_complement_inplace(ascii: &mut [u8]) {
    let len = ascii.len();
    for i in 0..len / 2 {
        let j = len - 1 - i;
        let c_i = COMPLEMENT_2BIT[ASCII_TO_2BIT[ascii[i] as usize] as usize];
        let c_j = COMPLEMENT_2BIT[ASCII_TO_2BIT[ascii[j] as usize] as usize];
        ascii[i] = TWO_BIT_TO_ASCII[c_j as usize];
        ascii[j] = TWO_BIT_TO_ASCII[c_i as usize];
    }
    if len % 2 == 1 {
        let mid = len / 2;
        let c_mid = COMPLEMENT_2BIT[ASCII_TO_2BIT[ascii[mid] as usize] as usize];
        ascii[mid] = TWO_BIT_TO_ASCII[c_mid as usize];
    }
}

/// Reverse-complements a packed 64-bit word of up to 32 bases.
/// Bases are assumed to be stored in bits `0..(2 * num_bases)` with base 0 at bits 0..1.
#[inline(always)]
pub fn rev_comp_word(word: u64, num_bases: usize) -> u64 {
    if num_bases == 0 {
        return 0;
    }
    debug_assert!(num_bases <= 32);

    let mut w = word.swap_bytes();
    w = ((w & 0x0f0f_0f0f_0f0f_0f0f) << 4) | ((w & 0xf0f0_f0f0_f0f0_f0f0) >> 4);
    w = ((w & 0x3333_3333_3333_3333) << 2) | ((w & 0xcccc_cccc_cccc_cccc) >> 2);
    w = !w;

    if num_bases < 32 {
        let shift = (32 - num_bases) * 2;
        let mask = (1u64 << (num_bases * 2)) - 1;
        (w >> shift) & mask
    } else {
        w
    }
}

/// Contiguous 2-bit sequence storage arena avoiding heap fragmentation.
///
/// Bases are tightly packed into `Vec<u64>` (32 bases per word).
/// Sequences are addressed by 0-based handle index `u32`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BitSeqArena {
    pub data: Vec<u64>,       // 32 bases per u64 word
    pub offsets: Vec<u64>,    // Base start index in `data`
    pub lengths: Vec<u32>,    // Sequence length in bp
}

impl BitSeqArena {
    /// Creates an empty sequence arena.
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            offsets: Vec::new(),
            lengths: Vec::new(),
        }
    }

    /// Creates an arena pre-allocated with capacity for `num_seqs` and `total_bases`.
    pub fn with_capacity(num_seqs: usize, total_bases: usize) -> Self {
        let num_words = (total_bases + 31) / 32;
        Self {
            data: Vec::with_capacity(num_words),
            offsets: Vec::with_capacity(num_seqs),
            lengths: Vec::with_capacity(num_seqs),
        }
    }

    /// Returns the number of sequences in the arena.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Returns true if the arena contains no sequences.
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Returns the total number of bases stored across all sequences.
    #[inline(always)]
    pub fn total_bases(&self) -> usize {
        if self.offsets.is_empty() {
            0
        } else {
            let last_idx = self.offsets.len() - 1;
            (self.offsets[last_idx] + self.lengths[last_idx] as u64) as usize
        }
    }

    /// Appends an ASCII DNA sequence into the arena. Returns the assigned sequence ID.
    pub fn push_seq(&mut self, ascii: &[u8]) -> u32 {
        let id = self.offsets.len() as u32;
        let start_pos = self.total_bases() as u64;
        let seq_len = ascii.len() as u32;

        self.offsets.push(start_pos);
        self.lengths.push(seq_len);

        if seq_len == 0 {
            return id;
        }

        let end_pos = start_pos + seq_len as u64;
        let required_words = ((end_pos + 31) / 32) as usize;
        if self.data.len() < required_words {
            self.data.resize(required_words, 0);
        }

        for (i, &b) in ascii.iter().enumerate() {
            let gpos = start_pos + i as u64;
            let word_idx = (gpos >> 5) as usize;
            let shift = (gpos & 31) << 1;
            let code = ASCII_TO_2BIT[b as usize] as u64;
            self.data[word_idx] |= code << shift;
        }

        id
    }

    /// Appends a raw 2-bit sequence slice (each byte in 0..3) into the arena.
    pub fn push_2bit(&mut self, bases: &[u8]) -> u32 {
        let id = self.offsets.len() as u32;
        let start_pos = self.total_bases() as u64;
        let seq_len = bases.len() as u32;

        self.offsets.push(start_pos);
        self.lengths.push(seq_len);

        if seq_len == 0 {
            return id;
        }

        let end_pos = start_pos + seq_len as u64;
        let required_words = ((end_pos + 31) / 32) as usize;
        if self.data.len() < required_words {
            self.data.resize(required_words, 0);
        }

        for (i, &b) in bases.iter().enumerate() {
            let gpos = start_pos + i as u64;
            let word_idx = (gpos >> 5) as usize;
            let shift = (gpos & 31) << 1;
            let code = (b & 0b11) as u64;
            self.data[word_idx] |= code << shift;
        }

        id
    }

    /// Returns the length of sequence `id` in base pairs.
    #[inline(always)]
    pub fn get_len(&self, id: u32) -> u32 {
        self.lengths[id as usize]
    }

    /// Retrieves the 2-bit encoded base (0..3) at position `pos` of sequence `id`.
    #[inline(always)]
    pub fn get_base(&self, id: u32, pos: usize) -> u8 {
        debug_assert!(pos < self.lengths[id as usize] as usize);
        let gpos = self.offsets[id as usize] + pos as u64;
        let word_idx = (gpos >> 5) as usize;
        let shift = (gpos & 31) << 1;
        ((self.data[word_idx] >> shift) & 0b11) as u8
    }

    /// Retrieves the ASCII character (`b'A'`, `b'C'`, `b'G'`, `b'T'`) at `pos` of `id`.
    #[inline(always)]
    pub fn get_base_ascii(&self, id: u32, pos: usize) -> u8 {
        TWO_BIT_TO_ASCII[self.get_base(id, pos) as usize]
    }

    /// Reconstructs the full sequence `id` as uppercase ASCII `Vec<u8>`.
    pub fn extract_ascii(&self, id: u32) -> Vec<u8> {
        let len = self.lengths[id as usize] as usize;
        let mut buf = Vec::with_capacity(len);
        self.extract_ascii_into(id, &mut buf);
        buf
    }

    /// Reconstructs sequence `id` into an existing `buf`, avoiding heap allocations.
    #[inline(always)]
    pub fn extract_ascii_into(&self, id: u32, buf: &mut Vec<u8>) {
        buf.clear();
        let len = self.lengths[id as usize] as usize;
        buf.reserve(len);
        let start_pos = self.offsets[id as usize];

        for i in 0..len {
            let gpos = start_pos + i as u64;
            let word_idx = (gpos >> 5) as usize;
            let shift = (gpos & 31) << 1;
            let code = ((self.data[word_idx] >> shift) & 0b11) as usize;
            buf.push(TWO_BIT_TO_ASCII[code]);
        }
    }

    /// Reconstructs a subslice of sequence `id` as uppercase ASCII `Vec<u8>`.
    pub fn extract_subseq_ascii(&self, id: u32, start: usize, len: usize) -> Vec<u8> {
        let seq_len = self.lengths[id as usize] as usize;
        assert!(start + len <= seq_len, "Subsequence out of range: start={}, len={}, seq_len={}", start, len, seq_len);
        let mut buf = Vec::with_capacity(len);
        let start_pos = self.offsets[id as usize] + start as u64;

        for i in 0..len {
            let gpos = start_pos + i as u64;
            let word_idx = (gpos >> 5) as usize;
            let shift = (gpos & 31) << 1;
            let code = ((self.data[word_idx] >> shift) & 0b11) as usize;
            buf.push(TWO_BIT_TO_ASCII[code]);
        }
        buf
    }

    /// Extracts the reverse complement of sequence `id` as uppercase ASCII `Vec<u8>`.
    pub fn extract_rc_ascii(&self, id: u32) -> Vec<u8> {
        let len = self.lengths[id as usize] as usize;
        let mut buf = Vec::with_capacity(len);
        let start_pos = self.offsets[id as usize];

        for i in (0..len).rev() {
            let gpos = start_pos + i as u64;
            let word_idx = (gpos >> 5) as usize;
            let shift = (gpos & 31) << 1;
            let code = ((self.data[word_idx] >> shift) & 0b11) as usize;
            let comp_code = COMPLEMENT_2BIT[code] as usize;
            buf.push(TWO_BIT_TO_ASCII[comp_code]);
        }
        buf
    }

    /// Computes the Hamming distance (number of mismatches) between two sequence slices
    /// in the arena, up to length `len`.
    pub fn mismatch_count(&self, id1: u32, start1: usize, id2: u32, start2: usize, len: usize) -> usize {
        let mut mismatches = 0;
        for i in 0..len {
            if self.get_base(id1, start1 + i) != self.get_base(id2, start2 + i) {
                mismatches += 1;
            }
        }
        mismatches
    }
}

/// Packed Read Pool separating merged single-end long reads from unmerged paired-end reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadBank {
    pub merged_reads: BitSeqArena,   // Overlapping PE reads merged into 180-280bp reads
    pub pe_forward: BitSeqArena,     // Unmerged R1 reads
    pub pe_reverse: BitSeqArena,     // Unmerged R2 reads
}

impl ReadBank {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(num_reads: usize, avg_len: usize) -> Self {
        Self {
            merged_reads: BitSeqArena::with_capacity(num_reads / 2, (num_reads / 2) * avg_len * 2),
            pe_forward: BitSeqArena::with_capacity(num_reads / 2, (num_reads / 2) * avg_len),
            pe_reverse: BitSeqArena::with_capacity(num_reads / 2, (num_reads / 2) * avg_len),
        }
    }

    pub fn total_reads(&self) -> usize {
        self.merged_reads.len() + self.pe_forward.len() + self.pe_reverse.len()
    }

    pub fn total_bases(&self) -> usize {
        self.merged_reads.total_bases() + self.pe_forward.total_bases() + self.pe_reverse.total_bases()
    }

    /// Iterates over all reads of length >= k without allocating a new Vec per read.
    pub fn for_each_read<F: FnMut(&[u8])>(&self, k: usize, mut f: F) {
        let mut buf = Vec::with_capacity(300);
        for id in 0..self.merged_reads.len() as u32 {
            if self.merged_reads.get_len(id) as usize >= k {
                self.merged_reads.extract_ascii_into(id, &mut buf);
                f(&buf);
            }
        }
        for id in 0..self.pe_forward.len() as u32 {
            if self.pe_forward.get_len(id) as usize >= k {
                self.pe_forward.extract_ascii_into(id, &mut buf);
                f(&buf);
            }
        }
        for id in 0..self.pe_reverse.len() as u32 {
            if self.pe_reverse.get_len(id) as usize >= k {
                self.pe_reverse.extract_ascii_into(id, &mut buf);
                f(&buf);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_encoding() {
        assert_eq!(Base::from_ascii(b'A'), Base::A);
        assert_eq!(Base::from_ascii(b'c'), Base::C);
        assert_eq!(Base::from_ascii(b'G'), Base::G);
        assert_eq!(Base::from_ascii(b't'), Base::T);
        assert_eq!(Base::from_ascii(b'N'), Base::A);

        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
        assert_eq!(Base::T.complement(), Base::A);
    }

    #[test]
    fn test_rev_comp_ascii() {
        let seq = b"ACGTACGT";
        assert_eq!(reverse_complement_ascii(seq), b"ACGTACGT");

        let seq2 = b"AAAACCCGGT";
        // Reverse: TGGCCCAAAA -> Complement: ACCGGGTTTT
        assert_eq!(reverse_complement_ascii(seq2), b"ACCGGGTTTT");

        let mut in_place = seq2.to_vec();
        reverse_complement_inplace(&mut in_place);
        assert_eq!(in_place, b"ACCGGGTTTT");
    }

    #[test]
    fn test_rev_comp_word() {
        // 3 bases: A (00), C (01), G (10) -> word = 00 01 10 = 0b000110 = 6
        // RC: C G T -> 01 10 11 = 0b011011 = 27
        let w = 0b00_01_10u64;
        let rc_w = rev_comp_word(w, 3);
        assert_eq!(rc_w, 0b01_10_11u64);

        // 4 bases: A C G T -> 00 01 10 11 = 0b00011011
        // RC: A C G T -> 00 01 10 11
        let w4 = 0b00_01_10_11u64;
        let rc_w4 = rev_comp_word(w4, 4);
        assert_eq!(rc_w4, 0b00_01_10_11u64);
    }

    #[test]
    fn test_bitseq_arena_packing_unpacking() {
        let mut arena = BitSeqArena::new();
        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"TGCATGCA";
        let seq3 = b"AAGGCCTT";

        let id1 = arena.push_seq(seq1);
        let id2 = arena.push_seq(seq2);
        let id3 = arena.push_seq(seq3);

        assert_eq!(id1, 0);
        assert_eq!(id2, 1);
        assert_eq!(id3, 2);

        assert_eq!(arena.extract_ascii(id1), seq1);
        assert_eq!(arena.extract_ascii(id2), seq2);
        assert_eq!(arena.extract_ascii(id3), seq3);

        assert_eq!(arena.extract_rc_ascii(id1), reverse_complement_ascii(seq1));
        assert_eq!(arena.extract_rc_ascii(id2), reverse_complement_ascii(seq2));
        assert_eq!(arena.extract_rc_ascii(id3), reverse_complement_ascii(seq3));
    }

    #[test]
    fn test_random_dna_lengths_1_to_500() {
        // Verification Gate item 1:
        // Packing and unpacking random DNA sequences of lengths 1..500 bp matches
        // original strings and their reverse complements.
        let mut arena = BitSeqArena::new();
        let alphabet = [b'A', b'C', b'G', b'T'];
        let mut rng_state: u64 = 0x123456789ABCDEF;

        let mut next_rand = || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            rng_state
        };

        for len in 1..=500 {
            let mut seq = Vec::with_capacity(len);
            for _ in 0..len {
                let idx = (next_rand() % 4) as usize;
                seq.push(alphabet[idx]);
            }

            let id = arena.push_seq(&seq);
            let unpacked = arena.extract_ascii(id);
            assert_eq!(unpacked, seq, "Failed on length {}", len);

            let rc_expected = reverse_complement_ascii(&seq);
            let rc_unpacked = arena.extract_rc_ascii(id);
            assert_eq!(rc_unpacked, rc_expected, "Failed RC on length {}", len);

            // Verify individual base access
            for pos in 0..len {
                assert_eq!(arena.get_base_ascii(id, pos), seq[pos]);
            }
        }
    }
}
