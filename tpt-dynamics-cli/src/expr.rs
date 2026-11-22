//! Row-expression language for `analyze --where`.
//!
//! Grammar, from lowest to highest precedence:
//!
//! ```text
//! expr    := and ( "or" and )*
//! and     := not ( "and" not )*
//! not     := "not" not | cmp
//! cmp     := add ( ( "<" | "<=" | ">" | ">=" | "==" | "!=" ) add )?
//! add     := mul ( ( "+" | "-" ) mul )*
//! mul     := unary ( ( "*" | "/" ) unary )*
//! unary   := "-" unary | primary
//! primary := NUMBER | COLUMN | "(" expr ")"
//! ```
//!
//! Comparisons do not chain: `a < b < c` is a syntax error. Numbers are decimal
//! literals (`2`, `0.5`, `.5`, `1e-3`). Column names are identifiers (letters, digits
//! and `_`, not starting with a digit). The keywords `and`, `or` and `not` are
//! lowercase only.
//!
//! Evaluation happens in the lazy-frame engine of `tpt-dynamics-data` (see [`lower`]),
//! so arithmetic is `f64` and works column-wise. Division by zero is not an error:
//! it follows IEEE 754, giving ±infinity for x / 0 and NaN for 0 / 0. A comparison
//! with NaN is false, so `<`, `<=`, `>`, `>=` and `==` drop such rows; `!=` keeps them
//! because NaN is not equal to anything (IEEE `!=` is true for NaN).

use tpt_dynamics_data::lazy::{col, Expr as Lazy};

/// Nesting limit for parentheses and `not`, to keep recursion bounded.
const MAX_DEPTH: usize = 128;

/// A parsed expression.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Node {
    /// A numeric literal.
    Number(f64),
    /// A column reference.
    Column(String),
    /// Arithmetic negation.
    Neg(Box<Node>),
    /// Logical negation.
    Not(Box<Node>),
    /// A binary operator applied to two operands.
    Binary(BinKind, Box<Node>, Box<Node>),
}

/// Binary operators of the expression language.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BinKind {
    Add,
    Sub,
    Mul,
    Div,
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
    And,
    Or,
}

/// A syntax error with the 1-based column where it was found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ParseError {
    /// What went wrong.
    pub(crate) message: String,
    /// 1-based column in the expression text.
    pub(crate) column: usize,
}

/// A name that the expression refers to but the table does not have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct UnknownColumn(pub(crate) String);

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(f64),
    Ident(String),
    Op(&'static str),
    LParen,
    RParen,
    And,
    Or,
    Not,
}

impl Tok {
    /// Text of the token as the user wrote it, for error messages.
    fn describe(&self) -> String {
        match self {
            Self::Num(v) => format!("number {v}"),
            Self::Ident(name) => format!("column '{name}'"),
            Self::Op(op) => format!("'{op}'"),
            Self::LParen => "'('".into(),
            Self::RParen => "')'".into(),
            Self::And => "'and'".into(),
            Self::Or => "'or'".into(),
            Self::Not => "'not'".into(),
        }
    }
}

/// Splits the text into tokens, each with its 1-based column.
fn lex(text: &str) -> Result<Vec<(Tok, usize)>, ParseError> {
    let bytes = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        let column = i + 1;
        let fail = |message: String| ParseError { message, column };
        if c.is_ascii_whitespace() {
            i += 1;
        } else if let Some(tok) = single_char_token(c) {
            out.push((tok, column));
            i += 1;
        } else if matches!(c, b'<' | b'>' | b'=' | b'!') {
            let (op, len) = lex_comparison(bytes, i).map_err(fail)?;
            out.push((Tok::Op(op), column));
            i += len;
        } else if c.is_ascii_digit()
            || (c == b'.' && bytes.get(i + 1).is_some_and(u8::is_ascii_digit))
        {
            let (value, end) = lex_number(text, i).map_err(fail)?;
            out.push((Tok::Num(value), column));
            i = end;
        } else if c.is_ascii_alphabetic() || c == b'_' {
            let end = word_end(bytes, i);
            out.push((word_token(&text[i..end]), column));
            i = end;
        } else {
            let ch = text[i..].chars().next().unwrap_or('?');
            return Err(fail(format!("unexpected character '{ch}'")));
        }
    }
    Ok(out)
}

/// The one-byte tokens: parentheses and the arithmetic operators.
fn single_char_token(c: u8) -> Option<Tok> {
    Some(match c {
        b'(' => Tok::LParen,
        b')' => Tok::RParen,
        b'+' => Tok::Op("+"),
        b'-' => Tok::Op("-"),
        b'*' => Tok::Op("*"),
        b'/' => Tok::Op("/"),
        _ => return None,
    })
}

/// The comparison operator starting at `i`, with its length in bytes.
fn lex_comparison(bytes: &[u8], i: usize) -> Result<(&'static str, usize), String> {
    let next = bytes.get(i + 1).copied();
    Ok(match (bytes[i], next) {
        (b'<', Some(b'=')) => ("<=", 2),
        (b'<', _) => ("<", 1),
        (b'>', Some(b'=')) => (">=", 2),
        (b'>', _) => (">", 1),
        (b'=', Some(b'=')) => ("==", 2),
        (b'!', Some(b'=')) => ("!=", 2),
        (b'=', _) => return Err("use '==' for equality".into()),
        _ => return Err("use 'not' for negation (and '!=' for inequality)".into()),
    })
}

/// The number starting at `start`: its value and the index just after it.
fn lex_number(text: &str, start: usize) -> Result<(f64, usize), String> {
    let bytes = text.as_bytes();
    let mut end = start;
    while end < bytes.len() && (bytes[end].is_ascii_digit() || bytes[end] == b'.') {
        end += 1;
    }
    if end < bytes.len() && (bytes[end] == b'e' || bytes[end] == b'E') {
        // An exponent is taken only when digits follow it.
        let mut j = end + 1;
        if j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
            j += 1;
        }
        if j < bytes.len() && bytes[j].is_ascii_digit() {
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            end = j;
        }
    }
    let lexeme = &text[start..end];
    let value = lexeme
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("malformed number '{lexeme}'"))?;
    Ok((value, end))
}

/// Index just after the identifier or keyword starting at `start`.
fn word_end(bytes: &[u8], start: usize) -> usize {
    let mut end = start;
    while end < bytes.len() && (bytes[end].is_ascii_alphanumeric() || bytes[end] == b'_') {
        end += 1;
    }
    end
}

/// A keyword token, or an identifier for any other word.
fn word_token(word: &str) -> Tok {
    match word {
        "and" => Tok::And,
        "or" => Tok::Or,
        "not" => Tok::Not,
        _ => Tok::Ident(word.to_string()),
    }
}

/// Recursive-descent parser over the token list.
struct Parser {
    toks: Vec<(Tok, usize)>,
    pos: usize,
    depth: usize,
    /// Column reported for "unexpected end of expression".
    end_column: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|(t, _)| t)
    }

    fn column(&self) -> usize {
        self.toks.get(self.pos).map_or(self.end_column, |(_, c)| *c)
    }

    fn error(&self, message: impl Into<String>) -> ParseError {
        ParseError { message: message.into(), column: self.column() }
    }

    fn bump(&mut self) -> Option<Tok> {
        let tok = self.toks.get(self.pos).map(|(t, _)| t.clone());
        self.pos += 1;
        tok
    }

    fn nest(&mut self) -> Result<(), ParseError> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err(self.error("expression is nested too deeply"));
        }
        Ok(())
    }

    fn or_expr(&mut self) -> Result<Node, ParseError> {
        let mut lhs = self.and_expr()?;
        while self.peek() == Some(&Tok::Or) {
            self.bump();
            let rhs = self.and_expr()?;
            lhs = Node::Binary(BinKind::Or, Box::new(lhs), Box::new(rhs));
        }
        Ok(lhs)
    }

    fn and_expr(&mut self) -> Result<Node, ParseError> {
        let mut lhs = self.not_expr()?;
        while self.peek() == Some(&Tok::And) {
            self.bump();
            let rhs = self.not_expr()?;
            lhs = Node::Binary(BinKind::And, Box::new(lhs), Box::new(rhs));
        }
        Ok(lhs)
    }

    fn not_expr(&mut self) -> Result<Node, ParseError> {
        if self.peek() == Some(&Tok::Not) {
            self.bump();
            self.nest()?;
            let inner = self.not_expr()?;
            self.depth -= 1;
            return Ok(Node::Not(Box::new(inner)));
        }
        self.cmp_expr()
    }

    fn cmp_expr(&mut self) -> Result<Node, ParseError> {
        let lhs = self.add_expr()?;
        let kind = match self.peek() {
            Some(Tok::Op("<")) => BinKind::Lt,
            Some(Tok::Op("<=")) => BinKind::Le,
            Some(Tok::Op(">")) => BinKind::Gt,
            Some(Tok::Op(">=")) => BinKind::Ge,
            Some(Tok::Op("==")) => BinKind::Eq,
            Some(Tok::Op("!=")) => BinKind::Ne,
            _ => return Ok(lhs),
        };
        self.bump();
        let rhs = self.add_expr()?;
        if matches!(self.peek(), Some(Tok::Op("<" | "<=" | ">" | ">=" | "==" | "!="))) {
            return Err(self.error("comparisons cannot be chained; use 'and' between them"));
        }
        Ok(Node::Binary(kind, Box::new(lhs), Box::new(rhs)))
    }

    fn add_expr(&mut self) -> Result<Node, ParseError> {
        let mut lhs = self.mul_expr()?;
        loop {
            let kind = match self.peek() {
                Some(Tok::Op("+")) => BinKind::Add,
                Some(Tok::Op("-")) => BinKind::Sub,
                _ => return Ok(lhs),
            };
            self.bump();
            let rhs = self.mul_expr()?;
            lhs = Node::Binary(kind, Box::new(lhs), Box::new(rhs));
        }
    }

    fn mul_expr(&mut self) -> Result<Node, ParseError> {
        let mut lhs = self.unary_expr()?;
        loop {
            let kind = match self.peek() {
                Some(Tok::Op("*")) => BinKind::Mul,
                Some(Tok::Op("/")) => BinKind::Div,
                _ => return Ok(lhs),
            };
            self.bump();
            let rhs = self.unary_expr()?;
            lhs = Node::Binary(kind, Box::new(lhs), Box::new(rhs));
        }
    }

    fn unary_expr(&mut self) -> Result<Node, ParseError> {
        if self.peek() == Some(&Tok::Op("-")) {
            self.bump();
            return Ok(Node::Neg(Box::new(self.unary_expr()?)));
        }
        self.primary()
    }

    fn primary(&mut self) -> Result<Node, ParseError> {
        match self.bump() {
            Some(Tok::Num(v)) => Ok(Node::Number(v)),
            Some(Tok::Ident(name)) => Ok(Node::Column(name)),
            Some(Tok::LParen) => {
                self.nest()?;
                let inner = self.or_expr()?;
                if self.peek() != Some(&Tok::RParen) {
                    return Err(self.error("missing ')'"));
                }
                self.bump();
                self.depth -= 1;
                Ok(inner)
            }
            Some(other) => {
                self.pos -= 1;
                Err(self.error(format!("expected a value, found {}", other.describe())))
            }
            None => Err(self.error("expression ends too early; expected a value")),
        }
    }
}

/// Parses an expression. The whole text must be consumed.
///
/// # Errors
/// Returns a [`ParseError`] with the column of the first problem.
pub(crate) fn parse(text: &str) -> Result<Node, ParseError> {
    let toks = lex(text)?;
    if toks.is_empty() {
        return Err(ParseError { message: "expression is empty".into(), column: 1 });
    }
    let end_column = text.len() + 1;
    let mut parser = Parser { toks, pos: 0, depth: 0, end_column };
    let node = parser.or_expr()?;
    if let Some(tok) = parser.peek() {
        let msg = format!("unexpected {} after the expression", tok.describe());
        return Err(parser.error(msg));
    }
    Ok(node)
}

/// Translates a parsed expression into the lazy engine's expression type.
///
/// `has_column` says whether the table has a column of that name. The first
/// unknown column found is reported.
///
/// # Errors
/// Returns [`UnknownColumn`] when the expression names a missing column.
pub(crate) fn lower(node: &Node, has_column: &dyn Fn(&str) -> bool) -> Result<Lazy, UnknownColumn> {
    Ok(match node {
        Node::Number(v) => Lazy::from(*v),
        Node::Column(name) => {
            if !has_column(name) {
                return Err(UnknownColumn(name.clone()));
            }
            col(name.as_str())
        }
        Node::Neg(a) => -lower(a, has_column)?,
        Node::Not(a) => !lower(a, has_column)?,
        Node::Binary(op, a, b) => {
            let l = lower(a, has_column)?;
            let r = lower(b, has_column)?;
            match op {
                BinKind::Add => l + r,
                BinKind::Sub => l - r,
                BinKind::Mul => l * r,
                BinKind::Div => l / r,
                BinKind::Lt => l.lt(r),
                BinKind::Le => l.le(r),
                BinKind::Gt => l.gt(r),
                BinKind::Ge => l.ge(r),
                BinKind::Eq => l.eq_to(r),
                BinKind::Ne => l.ne_to(r),
                BinKind::And => l.and(r),
                BinKind::Or => l.or(r),
            }
        }
    })
}

/// A value during row-wise evaluation (tests only; the engine evaluates column-wise).
#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Value {
    Num(f64),
    Bool(bool),
}

/// Evaluates one row. `lookup` returns a column's value in that row.
///
/// This mirrors the semantics the lazy engine implements, so the unit tests can check
/// precedence and NaN behaviour without a table.
#[cfg(test)]
#[allow(clippy::float_cmp)] // exact IEEE comparison is what the NaN tests check
pub(crate) fn eval(node: &Node, lookup: &dyn Fn(&str) -> Option<f64>) -> Result<Value, String> {
    let num = |v: Value| match v {
        Value::Num(x) => Ok(x),
        Value::Bool(_) => Err("expected a number".to_string()),
    };
    let boolean = |v: Value| match v {
        Value::Bool(b) => Ok(b),
        Value::Num(_) => Err("expected a boolean".to_string()),
    };
    Ok(match node {
        Node::Number(v) => Value::Num(*v),
        Node::Column(name) => Value::Num(lookup(name).ok_or_else(|| format!("no column {name}"))?),
        Node::Neg(a) => Value::Num(-num(eval(a, lookup)?)?),
        Node::Not(a) => Value::Bool(!boolean(eval(a, lookup)?)?),
        Node::Binary(op, a, b) => {
            let l = eval(a, lookup)?;
            let r = eval(b, lookup)?;
            match op {
                BinKind::Add => Value::Num(num(l)? + num(r)?),
                BinKind::Sub => Value::Num(num(l)? - num(r)?),
                BinKind::Mul => Value::Num(num(l)? * num(r)?),
                BinKind::Div => Value::Num(num(l)? / num(r)?),
                BinKind::Lt => Value::Bool(num(l)? < num(r)?),
                BinKind::Le => Value::Bool(num(l)? <= num(r)?),
                BinKind::Gt => Value::Bool(num(l)? > num(r)?),
                BinKind::Ge => Value::Bool(num(l)? >= num(r)?),
                BinKind::Eq => Value::Bool(num(l)? == num(r)?),
                BinKind::Ne => Value::Bool(num(l)? != num(r)?),
                BinKind::And => Value::Bool(boolean(l)? && boolean(r)?),
                BinKind::Or => Value::Bool(boolean(l)? || boolean(r)?),
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn value(text: &str) -> Value {
        let node = parse(text).unwrap_or_else(|e| panic!("{text}: {e:?}"));
        eval(&node, &|name| match name {
            "x" => Some(3.0),
            "y" => Some(-2.0),
            _ => None,
        })
        .unwrap_or_else(|e| panic!("{text}: {e}"))
    }

    fn num(text: &str) -> f64 {
        match value(text) {
            Value::Num(v) => v,
            Value::Bool(_) => panic!("{text} gave a boolean"),
        }
    }

    fn truth(text: &str) -> bool {
        match value(text) {
            Value::Bool(b) => b,
            Value::Num(_) => panic!("{text} gave a number"),
        }
    }

    #[test]
    fn multiplication_binds_tighter_than_addition() {
        assert_eq!(num("1 + 2 * 3"), 7.0);
        assert_eq!(num("8 - 4 / 2"), 6.0);
        assert_eq!(num("2 * 3 - 1"), 5.0);
    }

    #[test]
    fn parentheses_override_precedence() {
        assert_eq!(num("(1 + 2) * 3"), 9.0);
        assert_eq!(num("((x))*2"), 6.0);
        assert_eq!(num("-(x + 1)"), -4.0);
    }

    #[test]
    fn unary_minus_binds_to_the_operand() {
        assert_eq!(num("-x * 2"), -6.0);
        assert_eq!(num("2 - -x"), 5.0);
    }

    #[test]
    fn numeric_literal_forms() {
        assert_eq!(num("1e3"), 1000.0);
        assert_eq!(num(".5 + 1.5E-1"), 0.65);
        assert_eq!(num("2e+2"), 200.0);
    }

    #[test]
    fn and_binds_tighter_than_or() {
        // true or (false and false) = true; (true or false) and false = false.
        assert!(truth("1 < 2 or 3 < 2 and 3 < 2"));
        assert!(!truth("(1 < 2 or 3 < 2) and 3 < 2"));
    }

    #[test]
    fn not_applies_to_the_comparison_on_its_right() {
        assert!(truth("not 2 < 1"));
        assert!(!truth("not 1 < 2 and 3 == 3"));
        assert!(truth("not 1 == 2 and 3 == 3"));
    }

    #[test]
    fn comparisons_with_columns() {
        assert!(truth("x > 2 and y < 0"));
        assert!(!truth("x <= 2"));
        assert!(truth("x != y"));
        assert!(truth("x * 2 == 6"));
        assert!(truth("y >= -2"));
    }

    #[test]
    fn division_by_zero_gives_non_finite_values_that_compare_false() {
        let zero_over_zero = parse("0 / 0").unwrap();
        let lookup = |_: &str| None;
        let nan = match eval(&zero_over_zero, &lookup).unwrap() {
            Value::Num(v) => v,
            Value::Bool(_) => panic!("expected a number"),
        };
        assert!(nan.is_nan());
        assert!(!truth("0 / 0 > 1"));
        assert!(!truth("0 / 0 <= 1"));
        assert!(!truth("0 / 0 == 0 / 0"));
        assert!(truth("1 / 0 > 1000"));
    }

    #[test]
    fn syntax_errors_report_the_column() {
        let err = parse("1 +").unwrap_err();
        assert_eq!(err.column, 4);
        assert!(err.message.contains("ends too early"));

        let err = parse("pt > 2 extra").unwrap_err();
        assert_eq!(err.column, 8);
        assert!(err.message.contains("column 'extra'"), "{err:?}");

        let err = parse("(1 + 2").unwrap_err();
        assert!(err.message.contains("missing ')'"));

        let err = parse("x = 2").unwrap_err();
        assert!(err.message.contains("'=='"));
        assert_eq!(err.column, 3);

        assert!(parse("").unwrap_err().message.contains("empty"));
        assert!(parse("1 < 2 < 3").unwrap_err().message.contains("chained"));
        assert!(parse("1 @ 2").unwrap_err().message.contains("unexpected character '@'"));
        assert!(parse("1..2").is_err());
    }

    #[test]
    fn keywords_are_lowercase_only() {
        // "AND" is an identifier, so this is a column reference followed by a stray word.
        assert!(parse("1 < 2 AND 3 < 4").is_err());
    }

    #[test]
    fn lowering_rejects_unknown_columns() {
        let node = parse("pt > 20 and eta < 2").unwrap();
        let has = |name: &str| name == "pt";
        assert_eq!(lower(&node, &has).unwrap_err(), UnknownColumn("eta".into()));
        let node = parse("pt > 20").unwrap();
        assert!(lower(&node, &has).is_ok());
    }

    #[test]
    fn deep_nesting_is_rejected() {
        let text = format!("{}1{}", "(".repeat(200), ")".repeat(200));
        assert!(parse(&text).unwrap_err().message.contains("nested too deeply"));
    }
}
