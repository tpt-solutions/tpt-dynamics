//! CRC-32 (IEEE 802.3 polynomial, reflected, as used by zip and PNG).
//!
//! A small bitwise implementation kept in-crate so that checkpoint files do not
//! depend on an external checksum crate. Speed is not critical for checkpoints.

/// Reflected IEEE polynomial.
const POLY: u32 = 0xEDB8_8320;

/// Computes the CRC-32 of `data`.
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            let mask = (crc & 1).wrapping_neg();
            crc = (crc >> 1) ^ (POLY & mask);
        }
    }
    !crc
}

#[cfg(test)]
mod tests {
    use super::crc32;

    #[test]
    fn standard_check_value() {
        // Check value of the CRC-32/ISO-HDLC algorithm.
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc32(b""), 0);
    }
}
