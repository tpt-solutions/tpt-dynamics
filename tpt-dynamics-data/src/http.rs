//! Byte-range reads over plain HTTP/1.1, for [`crate::rowstore`].
//!
//! Each read is one `GET` with a `Range` header, answered with `206 Partial
//! Content`. The client is built on `std::net` only, so there is no TLS: the
//! URL must be `http://`. Serve data from a trusted network or put it behind a
//! TLS-terminating proxy you control.
//!
//! The client sends a generic `User-Agent` and no cookies or credentials.

use std::io::{Read, Write};
use std::net::TcpStream;
use std::time::Duration;

use crate::lazy::LazyError;
use crate::rowstore::RangeRead;

const TIMEOUT: Duration = Duration::from_secs(60);
/// Largest response header accepted, in bytes.
const MAX_HEADER: usize = 64 * 1024;

/// A parsed `http://host[:port]/path` URL.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Url {
    host: String,
    port: u16,
    path: String,
}

fn invalid(msg: impl Into<String>) -> LazyError {
    LazyError::Invalid(format!("http: {}", msg.into()))
}

fn is_plain_token(s: &str) -> bool {
    s.bytes().all(|b| b.is_ascii_graphic())
}

fn parse_url(url: &str) -> Result<Url, LazyError> {
    if url.starts_with("https://") {
        return Err(invalid("https is not supported (no TLS dependency); use http://"));
    }
    let rest = url.strip_prefix("http://").ok_or_else(|| invalid("URL must start with http://"))?;
    let (authority, path) = match rest.find('/') {
        Some(i) => (&rest[..i], &rest[i..]),
        None => (rest, "/"),
    };
    let (host, port) = match authority.rsplit_once(':') {
        Some((h, p)) if !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()) => {
            let port = p.parse::<u16>().map_err(|_| invalid(format!("bad port in '{url}'")))?;
            (h, port)
        }
        _ => (authority, 80),
    };
    if host.is_empty() || !is_plain_token(host) || !is_plain_token(path) {
        return Err(invalid(format!("malformed URL '{url}'")));
    }
    Ok(Url { host: host.to_string(), port, path: path.to_string() })
}

struct Response {
    status: u16,
    headers: Vec<(String, String)>,
    body: Vec<u8>,
}

impl Response {
    fn header(&self, name: &str) -> Option<&str> {
        self.headers.iter().find(|(n, _)| n == name).map(|(_, v)| v.as_str())
    }
}

/// One `GET` with a `Range` header. Returns the status, headers and body.
fn request(url: &Url, range: &str) -> Result<Response, LazyError> {
    let mut stream = TcpStream::connect((url.host.as_str(), url.port))
        .map_err(|e| invalid(format!("connect to {}:{}: {e}", url.host, url.port)))?;
    stream.set_read_timeout(Some(TIMEOUT)).map_err(|e| invalid(e.to_string()))?;
    stream.set_write_timeout(Some(TIMEOUT)).map_err(|e| invalid(e.to_string()))?;
    let host_header =
        if url.port == 80 { url.host.clone() } else { format!("{}:{}", url.host, url.port) };
    let request = format!(
        "GET {} HTTP/1.1\r\nHost: {host_header}\r\nRange: {range}\r\nConnection: close\r\nUser-Agent: tpt-dynamics\r\n\r\n",
        url.path
    );
    stream.write_all(request.as_bytes()).map_err(|e| invalid(format!("send: {e}")))?;

    // Read until the blank line that ends the header block.
    let mut buf = Vec::new();
    let split = loop {
        if let Some(pos) = buf.windows(4).position(|w| w == b"\r\n\r\n") {
            break pos;
        }
        if buf.len() > MAX_HEADER {
            return Err(invalid("response header too large"));
        }
        let mut chunk = [0u8; 4096];
        let n = stream.read(&mut chunk).map_err(|e| invalid(format!("receive: {e}")))?;
        if n == 0 {
            return Err(invalid("connection closed inside the response header"));
        }
        buf.extend_from_slice(&chunk[..n]);
    };
    let head = std::str::from_utf8(&buf[..split]).map_err(|_| invalid("header is not ASCII"))?;
    let mut lines = head.split("\r\n");
    let status_line = lines.next().unwrap_or("");
    let status = status_line
        .split_whitespace()
        .nth(1)
        .and_then(|s| s.parse::<u16>().ok())
        .ok_or_else(|| invalid(format!("bad status line '{status_line}'")))?;
    let headers: Vec<(String, String)> = lines
        .filter_map(|line| line.split_once(':'))
        .map(|(n, v)| (n.trim().to_ascii_lowercase(), v.trim().to_string()))
        .collect();
    let response_headers = Response { status, headers, body: Vec::new() };
    if response_headers
        .header("transfer-encoding")
        .is_some_and(|v| v.to_ascii_lowercase().contains("chunked"))
    {
        return Err(invalid("chunked transfer encoding is not supported for range reads"));
    }
    let leftover = buf.split_off(split + 4);
    let body = match response_headers.header("content-length").map(str::parse::<usize>) {
        Some(Ok(len)) => {
            let mut body = leftover;
            if body.len() > len {
                body.truncate(len);
            }
            let start = body.len();
            body.resize(len, 0);
            stream
                .read_exact(&mut body[start..])
                .map_err(|e| invalid(format!("body shorter than Content-Length: {e}")))?;
            body
        }
        Some(Err(_)) => return Err(invalid("bad Content-Length")),
        None => {
            // No length: the server must close the connection to end the body.
            let mut body = leftover;
            stream.read_to_end(&mut body).map_err(|e| invalid(format!("receive: {e}")))?;
            body
        }
    };
    Ok(Response { body, ..response_headers })
}

/// A remote byte sequence read with HTTP `Range` requests.
#[derive(Clone, Debug)]
pub struct HttpRange {
    url: Url,
    total: u64,
}

impl HttpRange {
    /// Opens the resource at `url` and learns its length with a one-byte range
    /// request.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the URL is not plain HTTP, the server
    /// does not answer range requests with `206`, or the length is missing.
    pub fn new(url: &str) -> Result<Self, LazyError> {
        let url = parse_url(url)?;
        let probe = request(&url, "bytes=0-0")?;
        if probe.status != 206 {
            return Err(invalid(format!(
                "server answered the range request with status {} (206 needed)",
                probe.status
            )));
        }
        let total = probe
            .header("content-range")
            .and_then(|v| v.rsplit_once('/'))
            .and_then(|(_, len)| len.trim().parse::<u64>().ok())
            .ok_or_else(|| invalid("Content-Range does not give the total length"))?;
        Ok(Self { url, total })
    }
}

impl RangeRead for HttpRange {
    fn byte_len(&self) -> Result<u64, LazyError> {
        Ok(self.total)
    }

    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), LazyError> {
        if buf.is_empty() {
            return Ok(());
        }
        let end = offset
            .checked_add(buf.len() as u64 - 1)
            .filter(|&e| e < self.total)
            .ok_or_else(|| invalid("read past the end of the resource"))?;
        let resp = request(&self.url, &format!("bytes={offset}-{end}"))?;
        if resp.status != 206 {
            return Err(invalid(format!("status {} for a range request", resp.status)));
        }
        if resp.body.len() != buf.len() {
            return Err(invalid(format!(
                "range returned {} bytes, {} requested",
                resp.body.len(),
                buf.len()
            )));
        }
        buf.copy_from_slice(&resp.body);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_host_port_and_path() {
        let u = parse_url("http://example.org:8080/data/events.rows").unwrap();
        assert_eq!(
            u,
            Url { host: "example.org".into(), port: 8080, path: "/data/events.rows".into() }
        );
        let u = parse_url("http://127.0.0.1").unwrap();
        assert_eq!(u, Url { host: "127.0.0.1".into(), port: 80, path: "/".into() });
    }

    #[test]
    fn rejects_https_and_malformed_urls() {
        assert!(parse_url("https://example.org/x").is_err());
        assert!(parse_url("ftp://example.org/x").is_err());
        assert!(parse_url("http:///x").is_err());
        assert!(parse_url("http://example.org/a b").is_err());
        assert!(parse_url("http://exa\r\nmple.org/").is_err());
        assert!(parse_url("http://example.org:99999/").is_err());
    }
}
