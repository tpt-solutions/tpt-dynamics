//! CRC-32 (IEEE 802.3, reflected polynomial `0xEDB88320`), implemented in-crate.
//!
//! This is the same checksum used by zlib, gzip and PNG. The lookup table is
//! built at compile time.

const POLY: u32 = 0xEDB8_8320;

const TABLE: [u32; 256] = build_table();

const fn build_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut crc = i as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = if crc & 1 != 0 { (crc >> 1) ^ POLY } else { crc >> 1 };
            bit += 1;
        }
        table[i] = crc;
        i += 1;
    }
    table
}

/// Incremental CRC-32 state.
///
/// ```
/// use tpt_dynamics_data::crc32::Crc32;
/// let mut c = Crc32::new();
/// c.update(b"1234");
/// c.update(b"56789");
/// assert_eq!(c.finish(), 0xCBF4_3926);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Crc32 {
    state: u32,
}

impl Default for Crc32 {
    fn default() -> Self {
        Self::new()
    }
}

impl Crc32 {
    /// Start a new checksum.
    #[inline]
    pub fn new() -> Self {
        Self { state: 0xFFFF_FFFF }
    }

    /// Feed bytes into the checksum.
    #[inline]
    pub fn update(&mut self, data: &[u8]) {
        let mut crc = self.state;
        for &b in data {
            crc = TABLE[((crc ^ u32::from(b)) & 0xFF) as usize] ^ (crc >> 8);
        }
        self.state = crc;
    }

    /// Final checksum value. The state is not consumed, so more data may be added
    /// afterwards.
    #[inline]
    pub fn finish(&self) -> u32 {
        !self.state
    }
}

/// CRC-32 of a byte slice in one call.
#[inline]
pub fn checksum(data: &[u8]) -> u32 {
    let mut c = Crc32::new();
    c.update(data);
    c.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_check_value() {
        assert_eq!(checksum(b"123456789"), 0xCBF4_3926);
        assert_eq!(checksum(b""), 0);
    }

    #[test]
    fn incremental_matches_oneshot() {
        let data: Vec<u8> = (0..=255u8).cycle().take(5000).collect();
        let mut c = Crc32::new();
        for chunk in data.chunks(77) {
            c.update(chunk);
        }
        assert_eq!(c.finish(), checksum(&data));
    }

    #[test]
    fn detects_single_bit_flip() {
        let mut data = b"detector block payload".to_vec();
        let good = checksum(&data);
        data[5] ^= 0x10;
        assert_ne!(checksum(&data), good);
    }
}
