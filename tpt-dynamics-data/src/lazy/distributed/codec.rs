//! Binary encoding of the values that cross the wire: expressions, aggregate
//! definitions, partial group states and data frames.
//!
//! Decoding reads untrusted bytes, so every length is checked against the bytes
//! that remain before anything is allocated, and expression nesting is capped.

use std::collections::HashMap;

use crate::lazy::exec::{Acc, KeyVal, Partial};
use crate::lazy::{
    AggExpr, AggKind, BinOp, Column, DType, DataFrame, Expr, Func, LazyError, Scalar, UnOp,
};

/// Maximum expression nesting accepted by the decoder.
const MAX_DEPTH: usize = 256;

const BINOPS: [BinOp; 12] = [
    BinOp::Add,
    BinOp::Sub,
    BinOp::Mul,
    BinOp::Div,
    BinOp::Eq,
    BinOp::Ne,
    BinOp::Lt,
    BinOp::Le,
    BinOp::Gt,
    BinOp::Ge,
    BinOp::And,
    BinOp::Or,
];
const UNOPS: [UnOp; 2] = [UnOp::Neg, UnOp::Not];
const FUNCS: [Func; 3] = [Func::InvariantMass, Func::DeltaR, Func::TransverseMomentum];
const AGG_KINDS: [AggKind; 5] =
    [AggKind::Count, AggKind::Sum, AggKind::Mean, AggKind::Min, AggKind::Max];

fn invalid(msg: impl Into<String>) -> LazyError {
    LazyError::Invalid(format!("wire: {}", msg.into()))
}

fn code_of<T: PartialEq>(table: &[T], value: &T) -> u8 {
    // Every value is in its table, so the position always exists.
    table.iter().position(|x| x == value).unwrap_or(0) as u8
}

fn from_code<T: Copy>(table: &[T], code: u8, what: &str) -> Result<T, LazyError> {
    table.get(code as usize).copied().ok_or_else(|| invalid(format!("unknown {what} code {code}")))
}

/// Append-only byte buffer.
#[derive(Default)]
pub(crate) struct Writer {
    buf: Vec<u8>,
}

impl Writer {
    pub(crate) fn finish(self) -> Vec<u8> {
        self.buf
    }

    pub(crate) fn u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    pub(crate) fn u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    pub(crate) fn u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    pub(crate) fn i64(&mut self, v: i64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    pub(crate) fn f64(&mut self, v: f64) {
        self.u64(v.to_bits());
    }

    pub(crate) fn len(&mut self, n: usize) {
        self.u32(u32::try_from(n).expect("wire: length exceeds u32"));
    }

    pub(crate) fn str(&mut self, s: &str) {
        self.len(s.len());
        self.buf.extend_from_slice(s.as_bytes());
    }
}

/// Bounds-checked cursor over received bytes.
pub(crate) struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    /// Fails unless every byte has been consumed.
    pub(crate) fn finish(&self) -> Result<(), LazyError> {
        if self.pos == self.data.len() {
            Ok(())
        } else {
            Err(invalid("trailing bytes after message"))
        }
    }

    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8], LazyError> {
        let end = self.pos.checked_add(n).filter(|&e| e <= self.data.len());
        let end = end.ok_or_else(|| invalid("message truncated"))?;
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    pub(crate) fn u8(&mut self) -> Result<u8, LazyError> {
        Ok(self.take(1)?[0])
    }

    pub(crate) fn u32(&mut self) -> Result<u32, LazyError> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes")))
    }

    pub(crate) fn u64(&mut self) -> Result<u64, LazyError> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().expect("8 bytes")))
    }

    pub(crate) fn i64(&mut self) -> Result<i64, LazyError> {
        Ok(i64::from_le_bytes(self.take(8)?.try_into().expect("8 bytes")))
    }

    pub(crate) fn f64(&mut self) -> Result<f64, LazyError> {
        Ok(f64::from_bits(self.u64()?))
    }

    /// Reads an element count, rejecting counts that cannot fit in what remains
    /// (each element takes at least `min_bytes`).
    pub(crate) fn count(&mut self, min_bytes: usize) -> Result<usize, LazyError> {
        let n = self.u32()? as usize;
        if n.saturating_mul(min_bytes) > self.data.len() - self.pos {
            return Err(invalid("element count exceeds message size"));
        }
        Ok(n)
    }

    pub(crate) fn str(&mut self) -> Result<String, LazyError> {
        let n = self.count(1)?;
        String::from_utf8(self.take(n)?.to_vec()).map_err(|_| invalid("string is not UTF-8"))
    }
}

fn put_scalar(w: &mut Writer, s: &Scalar) {
    match s {
        Scalar::F64(v) => {
            w.u8(0);
            w.f64(*v);
        }
        Scalar::I64(v) => {
            w.u8(1);
            w.i64(*v);
        }
        Scalar::Bool(v) => {
            w.u8(2);
            w.u8(u8::from(*v));
        }
    }
}

fn get_scalar(r: &mut Reader<'_>) -> Result<Scalar, LazyError> {
    Ok(match r.u8()? {
        0 => Scalar::F64(r.f64()?),
        1 => Scalar::I64(r.i64()?),
        2 => Scalar::Bool(r.u8()? != 0),
        other => return Err(invalid(format!("unknown scalar tag {other}"))),
    })
}

pub(crate) fn put_expr(w: &mut Writer, e: &Expr) {
    match e {
        Expr::Col(name) => {
            w.u8(0);
            w.str(name);
        }
        Expr::Lit(s) => {
            w.u8(1);
            put_scalar(w, s);
        }
        Expr::Binary { op, lhs, rhs } => {
            w.u8(2);
            w.u8(code_of(&BINOPS, op));
            put_expr(w, lhs);
            put_expr(w, rhs);
        }
        Expr::Unary { op, arg } => {
            w.u8(3);
            w.u8(code_of(&UNOPS, op));
            put_expr(w, arg);
        }
        Expr::Call { func, args } => {
            w.u8(4);
            w.u8(code_of(&FUNCS, func));
            w.len(args.len());
            for a in args {
                put_expr(w, a);
            }
        }
        Expr::Alias { expr, name } => {
            w.u8(5);
            put_expr(w, expr);
            w.str(name);
        }
    }
}

pub(crate) fn get_expr(r: &mut Reader<'_>) -> Result<Expr, LazyError> {
    get_expr_depth(r, 0)
}

fn get_expr_depth(r: &mut Reader<'_>, depth: usize) -> Result<Expr, LazyError> {
    if depth > MAX_DEPTH {
        return Err(invalid("expression nested too deeply"));
    }
    let next = |r: &mut Reader<'_>| get_expr_depth(r, depth + 1);
    Ok(match r.u8()? {
        0 => Expr::Col(r.str()?),
        1 => Expr::Lit(get_scalar(r)?),
        2 => {
            let op = from_code(&BINOPS, r.u8()?, "binary operator")?;
            Expr::Binary { op, lhs: Box::new(next(r)?), rhs: Box::new(next(r)?) }
        }
        3 => {
            let op = from_code(&UNOPS, r.u8()?, "unary operator")?;
            Expr::Unary { op, arg: Box::new(next(r)?) }
        }
        4 => {
            let func = from_code(&FUNCS, r.u8()?, "function")?;
            let n = r.count(1)?;
            let args = (0..n).map(|_| next(r)).collect::<Result<Vec<_>, _>>()?;
            Expr::Call { func, args }
        }
        5 => {
            let expr = Box::new(next(r)?);
            Expr::Alias { expr, name: r.str()? }
        }
        other => return Err(invalid(format!("unknown expression tag {other}"))),
    })
}

pub(crate) fn put_agg(w: &mut Writer, a: &AggExpr) {
    w.u8(code_of(&AGG_KINDS, &a.kind()));
    match a.input() {
        Some(e) => {
            w.u8(1);
            put_expr(w, e);
        }
        None => w.u8(0),
    }
    w.str(a.name());
}

pub(crate) fn get_agg(r: &mut Reader<'_>) -> Result<AggExpr, LazyError> {
    let kind = from_code(&AGG_KINDS, r.u8()?, "aggregate")?;
    let input = match r.u8()? {
        0 => None,
        1 => Some(get_expr(r)?),
        other => return Err(invalid(format!("bad aggregate input flag {other}"))),
    };
    let name = r.str()?;
    Ok(AggExpr::from_parts(kind, input, name))
}

fn put_key(w: &mut Writer, k: &KeyVal) {
    match k {
        KeyVal::I64(v) => {
            w.u8(0);
            w.i64(*v);
        }
        KeyVal::Bool(v) => {
            w.u8(1);
            w.u8(u8::from(*v));
        }
    }
}

fn get_key(r: &mut Reader<'_>) -> Result<KeyVal, LazyError> {
    Ok(match r.u8()? {
        0 => KeyVal::I64(r.i64()?),
        1 => KeyVal::Bool(r.u8()? != 0),
        other => return Err(invalid(format!("unknown key tag {other}"))),
    })
}

fn put_acc(w: &mut Writer, a: &Acc) {
    match *a {
        Acc::Count(n) => {
            w.u8(0);
            w.i64(n);
        }
        Acc::SumI(n) => {
            w.u8(1);
            w.i64(n);
        }
        Acc::SumF(x) => {
            w.u8(2);
            w.f64(x);
        }
        Acc::Mean { sum, n } => {
            w.u8(3);
            w.f64(sum);
            w.i64(n);
        }
        Acc::MinF(x) => {
            w.u8(4);
            w.f64(x);
        }
        Acc::MaxF(x) => {
            w.u8(5);
            w.f64(x);
        }
        Acc::MinI(n) => {
            w.u8(6);
            w.i64(n);
        }
        Acc::MaxI(n) => {
            w.u8(7);
            w.i64(n);
        }
    }
}

fn get_acc(r: &mut Reader<'_>) -> Result<Acc, LazyError> {
    Ok(match r.u8()? {
        0 => Acc::Count(r.i64()?),
        1 => Acc::SumI(r.i64()?),
        2 => Acc::SumF(r.f64()?),
        3 => Acc::Mean { sum: r.f64()?, n: r.i64()? },
        4 => Acc::MinF(r.f64()?),
        5 => Acc::MaxF(r.f64()?),
        6 => Acc::MinI(r.i64()?),
        7 => Acc::MaxI(r.i64()?),
        other => return Err(invalid(format!("unknown accumulator tag {other}"))),
    })
}

/// Encodes a partial aggregate. Groups are written in key order so the bytes
/// do not depend on hash-map iteration order.
pub(crate) fn put_partial(w: &mut Writer, p: &Partial) {
    let mut groups: Vec<(&Vec<KeyVal>, &Vec<Acc>)> = p.iter().collect();
    groups.sort_by(|a, b| a.0.cmp(b.0));
    w.len(groups.len());
    for (key, accs) in groups {
        w.len(key.len());
        for k in key {
            put_key(w, k);
        }
        w.len(accs.len());
        for a in accs {
            put_acc(w, a);
        }
    }
}

pub(crate) fn get_partial(r: &mut Reader<'_>) -> Result<Partial, LazyError> {
    let groups = r.count(8)?;
    let mut map = HashMap::with_capacity(groups);
    for _ in 0..groups {
        let nkeys = r.count(2)?;
        let key = (0..nkeys).map(|_| get_key(r)).collect::<Result<Vec<_>, _>>()?;
        let naccs = r.count(9)?;
        let accs = (0..naccs).map(|_| get_acc(r)).collect::<Result<Vec<_>, _>>()?;
        if map.insert(key, accs).is_some() {
            return Err(invalid("duplicate group in partial aggregate"));
        }
    }
    Ok(map)
}

fn dtype_code(d: DType) -> u8 {
    match d {
        DType::I64 => 0,
        DType::F64 => 1,
        DType::Bool => 2,
    }
}

/// Encodes a frame: row count, then each named column with its type.
pub(crate) fn put_frame(w: &mut Writer, df: &DataFrame) -> Result<(), LazyError> {
    w.u64(df.len() as u64);
    w.len(df.names().len());
    for name in df.names() {
        let column = df.column(name).ok_or_else(|| LazyError::UnknownColumn(name.clone()))?;
        w.str(name);
        w.u8(dtype_code(column.dtype()));
        match column {
            Column::I64(v) => v.iter().for_each(|&x| w.i64(x)),
            Column::F64(v) => v.iter().for_each(|&x| w.f64(x)),
            Column::Bool(v) => v.iter().for_each(|&x| w.u8(u8::from(x))),
        }
    }
    Ok(())
}

pub(crate) fn get_frame(r: &mut Reader<'_>) -> Result<DataFrame, LazyError> {
    let rows = usize::try_from(r.u64()?).map_err(|_| invalid("row count too large"))?;
    let ncols = r.count(6)?;
    let mut columns = Vec::with_capacity(ncols);
    for _ in 0..ncols {
        let name = r.str()?;
        let column = match r.u8()? {
            0 => {
                let bytes = rows.checked_mul(8).ok_or_else(|| invalid("column too large"))?;
                let raw = r.take(bytes)?;
                Column::I64(raw.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)).collect())
            }
            1 => {
                let bytes = rows.checked_mul(8).ok_or_else(|| invalid("column too large"))?;
                let raw = r.take(bytes)?;
                Column::F64(
                    raw.as_chunks::<8>()
                        .0
                        .iter()
                        .map(|c| f64::from_bits(u64::from_le_bytes(*c)))
                        .collect(),
                )
            }
            2 => Column::Bool(r.take(rows)?.iter().map(|&b| b != 0).collect()),
            other => return Err(invalid(format!("unknown column type {other}"))),
        };
        columns.push((name, column));
    }
    DataFrame::with_rows(rows, columns)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::lazy::{col, invariant_mass, lit};
    use crate::rng::SplitMix64;

    fn encode_expr(e: &Expr) -> Vec<u8> {
        let mut w = Writer::default();
        put_expr(&mut w, e);
        w.finish()
    }

    fn sample_exprs() -> Vec<Expr> {
        vec![
            col("pt"),
            lit(1.5),
            lit(-3_i64),
            lit(true),
            -(col("x") * lit(2.0) + col("y")).alias("z"),
            !col("flag").and(col("g").gt(0_i64)).or(col("h").le(1.0)),
            invariant_mass(
                col("a"),
                col("b"),
                col("c"),
                col("d"),
                lit(1.0),
                lit(0.5),
                lit(0.25),
                lit(0.0),
            ),
        ]
    }

    #[test]
    fn expressions_round_trip() {
        for e in sample_exprs() {
            let bytes = encode_expr(&e);
            let mut r = Reader::new(&bytes);
            let back = get_expr(&mut r).unwrap();
            r.finish().unwrap();
            assert_eq!(back, e);
        }
    }

    #[test]
    fn partial_round_trips_every_accumulator() {
        let mut p: Partial = HashMap::new();
        p.insert(
            vec![KeyVal::I64(-1), KeyVal::Bool(true)],
            vec![
                Acc::Count(7),
                Acc::SumI(-9),
                Acc::SumF(2.5),
                Acc::Mean { sum: 10.0, n: 4 },
                Acc::MinF(-0.5),
                Acc::MaxF(f64::INFINITY),
                Acc::MinI(-3),
                Acc::MaxI(8),
            ],
        );
        let mut w = Writer::default();
        put_partial(&mut w, &p);
        let bytes = w.finish();
        let mut r = Reader::new(&bytes);
        let back = get_partial(&mut r).unwrap();
        r.finish().unwrap();
        assert_eq!(format!("{back:?}"), format!("{p:?}"));
    }

    #[test]
    fn frame_round_trips_every_column_type() {
        let df = DataFrame::new(vec![
            ("i", Column::I64(vec![1, -2, 3])),
            ("f", Column::F64(vec![0.1, f64::NAN, -0.0])),
            ("b", Column::Bool(vec![true, false, true])),
        ])
        .unwrap();
        let mut w = Writer::default();
        put_frame(&mut w, &df).unwrap();
        let bytes = w.finish();
        let mut r = Reader::new(&bytes);
        let back = get_frame(&mut r).unwrap();
        r.finish().unwrap();
        assert_eq!(back.names(), df.names());
        assert_eq!(back.column("i"), df.column("i"));
        assert_eq!(back.column("b"), df.column("b"));
        // NaN != NaN, so compare the bit patterns.
        let (Some(Column::F64(a)), Some(Column::F64(b))) = (back.column("f"), df.column("f"))
        else {
            panic!("f must be f64");
        };
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(a), bits(b));
    }

    #[test]
    fn deeply_nested_expression_is_rejected() {
        let mut w = Writer::default();
        for _ in 0..(MAX_DEPTH + 10) {
            w.u8(3);
            w.u8(0);
        }
        w.u8(0);
        w.str("x");
        let bytes = w.finish();
        assert!(get_expr(&mut Reader::new(&bytes)).is_err());
    }

    #[test]
    fn truncated_messages_are_errors_not_panics() {
        for e in sample_exprs() {
            let bytes = encode_expr(&e);
            for cut in 0..bytes.len() {
                let mut r = Reader::new(&bytes[..cut]);
                let res = get_expr(&mut r).and_then(|_| r.finish());
                assert!(res.is_err(), "prefix of length {cut} decoded");
            }
        }
    }

    #[test]
    fn random_bytes_never_panic() {
        let mut rng = SplitMix64::new(0xC0DE);
        for _ in 0..5_000 {
            let len = rng.below(96) as usize;
            let bytes: Vec<u8> = (0..len).map(|_| rng.next_u64() as u8).collect();
            let _ = get_expr(&mut Reader::new(&bytes));
            let _ = get_partial(&mut Reader::new(&bytes));
            let _ = get_frame(&mut Reader::new(&bytes));
            let _ = get_agg(&mut Reader::new(&bytes));
        }
    }
}
