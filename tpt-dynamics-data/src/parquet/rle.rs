//! Parquet RLE / bit-packed hybrid decoder (used for levels and dictionary
//! indices).

use super::ParquetError;

/// Reads a little-endian unsigned varint from `data` at `*pos`.
fn read_varint(data: &[u8], pos: &mut usize) -> Result<u64, ParquetError> {
    let mut out = 0u64;
    for shift in (0..64).step_by(7) {
        let b = *data
            .get(*pos)
            .ok_or_else(|| ParquetError::Format("truncated RLE run header".into()))?;
        *pos += 1;
        out |= u64::from(b & 0x7f) << shift;
        if b & 0x80 == 0 {
            return Ok(out);
        }
    }
    Err(ParquetError::Format("RLE run header too long".into()))
}

/// Decodes `count` values of `bit_width` bits each from the hybrid encoding.
pub(crate) fn decode_hybrid(
    data: &[u8],
    bit_width: u32,
    count: usize,
) -> Result<Vec<u32>, ParquetError> {
    if bit_width > 32 {
        return Err(ParquetError::Format(format!("RLE bit width {bit_width} exceeds 32")));
    }
    let mut out: Vec<u32> = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| ParquetError::Format(format!("cannot allocate {count} RLE values")))?;
    let mut pos = 0usize;
    let value_bytes = bit_width.div_ceil(8) as usize;
    while out.len() < count {
        let header = read_varint(data, &mut pos)?;
        if header & 1 == 1 {
            // Bit-packed run: groups of 8 values.
            let groups = (header >> 1) as usize;
            let nbytes = groups
                .checked_mul(bit_width as usize)
                .ok_or_else(|| ParquetError::Format("bit-packed run too large".into()))?;
            let end = pos
                .checked_add(nbytes)
                .filter(|&e| e <= data.len())
                .ok_or_else(|| ParquetError::Format("truncated bit-packed run".into()))?;
            let packed = &data[pos..end];
            pos = end;
            let take = groups.saturating_mul(8).min(count - out.len());
            for i in 0..take {
                out.push(unpack(packed, bit_width, i));
            }
        } else {
            // RLE run: one value repeated.
            let run = (header >> 1) as usize;
            let end = pos
                .checked_add(value_bytes)
                .filter(|&e| e <= data.len())
                .ok_or_else(|| ParquetError::Format("truncated RLE run value".into()))?;
            let mut v = 0u32;
            for (i, &b) in data[pos..end].iter().enumerate() {
                v |= u32::from(b) << (8 * i);
            }
            pos = end;
            let take = run.min(count - out.len());
            out.extend(std::iter::repeat_n(v, take));
        }
    }
    Ok(out)
}

/// Reads the `i`-th `bit_width`-bit value from an LSB-first packed buffer.
fn unpack(packed: &[u8], bit_width: u32, i: usize) -> u32 {
    let mut v = 0u32;
    for b in 0..bit_width as usize {
        let bit = i * bit_width as usize + b;
        if (packed[bit / 8] >> (bit % 8)) & 1 == 1 {
            v |= 1 << b;
        }
    }
    v
}

/// Decodes `count` PLAIN-encoded booleans (bit-packed, LSB first).
pub(crate) fn decode_bools(data: &[u8], count: usize) -> Result<Vec<bool>, ParquetError> {
    if data.len().saturating_mul(8) < count {
        return Err(ParquetError::Format("truncated boolean page".into()));
    }
    Ok((0..count).map(|i| unpack(data, 1, i) == 1).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rle_run_repeats_value() {
        // header = 5 << 1 (RLE, count 5); value 3 in one byte at bit width 3.
        assert_eq!(decode_hybrid(&[10, 3], 3, 5).unwrap(), vec![3; 5]);
    }

    #[test]
    fn bit_packed_run_is_lsb_first() {
        // header = (1 group << 1) | 1; 8 values of width 3: 0,1,2,...,7.
        let packed: [u8; 3] = [0b1000_1000, 0b1100_0110, 0b1111_1010];
        let mut data = vec![3u8];
        data.extend_from_slice(&packed);
        assert_eq!(decode_hybrid(&data, 3, 8).unwrap(), vec![0, 1, 2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn truncated_run_is_an_error() {
        assert!(decode_hybrid(&[10], 3, 5).is_err());
    }

    #[test]
    fn booleans_are_bit_packed() {
        assert_eq!(decode_bools(&[0b0000_0101], 3).unwrap(), vec![true, false, true]);
    }
}
