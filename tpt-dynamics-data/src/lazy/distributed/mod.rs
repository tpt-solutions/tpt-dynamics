//! Distributed execution of lazy plans across processes (todo Phase 7).
//!
//! A coordinator splits a source's rows into one contiguous range per worker.
//! Each worker reads its range, runs the row-wise stages and either returns the
//! filtered rows or, for a group-by, its unfinalised partial aggregate. The
//! coordinator merges the replies in worker order and finalises them, so the
//! result matches a local run of the same plan.
//!
//! Workers are plain TCP servers speaking a small length-prefixed binary
//! protocol (see `codec`). They resolve source names themselves, so the
//! source data must be reachable from every worker, for example on shared
//! storage such as a [`crate::rowstore::RowStore`] file.
//!
//! The protocol carries no authentication or encryption. Run workers only on a
//! trusted network.
//!
//! ```no_run
//! use std::net::TcpListener;
//! use std::sync::Arc;
//! use tpt_dynamics_data::lazy::distributed::{collect, serve, DistributedPlan};
//! use tpt_dynamics_data::lazy::{col, AggExpr};
//! use tpt_dynamics_data::rowstore::RowStore;
//!
//! // On each worker machine:
//! // let listener = TcpListener::bind("0.0.0.0:7000").unwrap();
//! // serve(&listener, &Arc::new(|_name| Ok(Arc::new(RowStore::open("events.rows")?) as _)), None).unwrap();
//!
//! // On the coordinator:
//! let store = RowStore::open("events.rows").unwrap();
//! let plan = DistributedPlan::new()
//!     .filter(col("pt").gt(20.0))
//!     .group_by_agg(["charge"], [AggExpr::count(), AggExpr::sum(col("pt")).alias("pt_sum")]);
//! let workers = ["10.0.0.1:7000".parse().unwrap(), "10.0.0.2:7000".parse().unwrap()];
//! let out = collect(&plan, "events", &store, &workers).unwrap();
//! # let _ = (out, Arc::new(0));
//! ```

pub(crate) mod codec;

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::ops::Range;
use std::sync::Arc;
use std::time::Duration;

use super::exec::{self, CollectOptions, Partial};
use super::expr::Expr;
use super::frame::{Column, DataFrame, LazyError, LazyFrame, Source};
use super::plan::{validate, AggExpr};
use codec::{Reader, Writer};

/// Protocol version. Peers with a different version refuse each other.
const PROTOCOL_VERSION: u8 = 1;
/// Largest accepted message, in bytes.
const MAX_MESSAGE: usize = 1 << 31;
/// Timeout for opening a connection to a worker.
const CONNECT_TIMEOUT: Duration = Duration::from_secs(10);

/// Maps a source name from the coordinator to a source, on a worker.
pub type ResolveFn = dyn Fn(&str) -> Result<Arc<dyn Source>, LazyError> + Send + Sync;
/// Shared handle to a [`ResolveFn`], as [`serve`] takes it.
pub type Resolver = Arc<ResolveFn>;

const MSG_JOB: u8 = 1;
const REPLY_FRAME: u8 = 0;
const REPLY_PARTIAL: u8 = 1;
const REPLY_ERROR: u8 = 2;

/// A row-wise operator in a [`DistributedPlan`].
#[derive(Clone, Debug, PartialEq)]
pub enum Stage {
    /// Keep the rows where the boolean expression holds.
    Filter(Expr),
    /// Replace the columns with the named expressions.
    Project(Vec<Expr>),
    /// Add or replace one named column.
    WithColumn(Expr),
}

/// A plan that can be shipped to workers: row-wise stages, then an optional
/// group-by aggregate. Sources are not part of the plan; the coordinator names
/// one and each worker resolves it.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct DistributedPlan {
    stages: Vec<Stage>,
    aggregate: Option<(Vec<String>, Vec<AggExpr>)>,
}

impl DistributedPlan {
    /// An empty plan that returns its source unchanged.
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a filter.
    #[must_use]
    pub fn filter(mut self, predicate: Expr) -> Self {
        self.stages.push(Stage::Filter(predicate));
        self
    }

    /// Appends a projection onto the named expressions.
    #[must_use]
    pub fn select(mut self, exprs: impl IntoIterator<Item = Expr>) -> Self {
        self.stages.push(Stage::Project(exprs.into_iter().collect()));
        self
    }

    /// Appends a with-column stage.
    #[must_use]
    pub fn with_column(mut self, expr: Expr) -> Self {
        self.stages.push(Stage::WithColumn(expr));
        self
    }

    /// Ends the plan with a group-by aggregate. Replaces any earlier aggregate.
    #[must_use]
    pub fn group_by_agg(
        mut self,
        keys: impl IntoIterator<Item = impl Into<String>>,
        aggs: impl IntoIterator<Item = AggExpr>,
    ) -> Self {
        self.aggregate =
            Some((keys.into_iter().map(Into::into).collect(), aggs.into_iter().collect()));
        self
    }

    /// The same plan as a local [`LazyFrame`] over `source`. Useful as a
    /// reference when checking a distributed result.
    pub fn to_lazy(&self, source: Arc<dyn Source>) -> LazyFrame {
        let mut lf = LazyFrame::scan_arc(source);
        for stage in &self.stages {
            lf = match stage {
                Stage::Filter(e) => lf.filter(e.clone()),
                Stage::Project(es) => lf.select(es.clone()),
                Stage::WithColumn(e) => lf.with_column(e.clone()),
            };
        }
        if let Some((keys, aggs)) = &self.aggregate {
            lf = lf.group_by(keys.clone()).agg(aggs.clone());
        }
        lf
    }
}

/// Runs `plan` over `source` on `workers` and returns the combined result.
///
/// `source_name` is the name each worker uses to look up the source, and
/// `source` supplies the schema and row count the coordinator needs.
///
/// # Errors
///
/// Returns [`LazyError::Invalid`] if there are no workers, a worker cannot be
/// reached, or a worker reports an error. Returns the usual plan validation
/// errors for an invalid plan.
pub fn collect(
    plan: &DistributedPlan,
    source_name: &str,
    source: &dyn Source,
    workers: &[SocketAddr],
) -> Result<DataFrame, LazyError> {
    if workers.is_empty() {
        return Err(LazyError::Invalid("no workers given".into()));
    }
    // The output schema depends only on the input schema, so validate the plan
    // against an empty frame with the same columns.
    let schema = source.schema();
    let empty = DataFrame::with_rows(
        0,
        schema.iter().map(|f| (f.name.clone(), Column::empty(f.dtype))).collect::<Vec<_>>(),
    )?;
    let local = plan.to_lazy(Arc::new(empty)).plan().clone();
    let out_schema = validate(&local)?;

    let ranges = split_rows(source.num_rows(), workers.len());
    let replies = std::thread::scope(|scope| {
        let handles: Vec<_> = workers
            .iter()
            .zip(ranges)
            .map(|(addr, rows)| {
                let job = encode_job(source_name, rows, plan);
                scope.spawn(move || send_job(*addr, &job))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .unwrap_or_else(|_| Err(LazyError::Invalid("worker thread panicked".into())))
            })
            .collect::<Result<Vec<_>, _>>()
    })?;

    // Replies are in worker order, so the merge order is fixed.
    if plan.aggregate.is_some() {
        let mut groups = Partial::new();
        for reply in replies {
            match reply {
                Reply::Partial(p) => exec::merge(&mut groups, p),
                Reply::Frame(_) => {
                    return Err(LazyError::Invalid("worker sent rows for an aggregate plan".into()))
                }
            }
        }
        exec::finalize(groups, &out_schema)
    } else {
        let mut parts = Vec::with_capacity(replies.len());
        for reply in replies {
            match reply {
                Reply::Frame(df) => parts.push(df),
                Reply::Partial(_) => {
                    return Err(LazyError::Invalid("worker sent a partial for a row plan".into()))
                }
            }
        }
        exec::concat(&out_schema, &parts)
    }
}

/// Serves jobs on `listener`, one connection at a time. Stops after `max_jobs`
/// connections, or runs until an accept error if `max_jobs` is `None`.
///
/// `resolve` maps a source name from the coordinator to a source. A job that
/// fails is answered with an error message and the server keeps going.
///
/// # Errors
///
/// Returns [`LazyError::Invalid`] if the listener cannot accept connections.
pub fn serve(
    listener: &TcpListener,
    resolve: &Resolver,
    max_jobs: Option<usize>,
) -> Result<(), LazyError> {
    let mut served = 0usize;
    loop {
        if max_jobs.is_some_and(|m| served >= m) {
            return Ok(());
        }
        let (mut stream, _) =
            listener.accept().map_err(|e| LazyError::Invalid(format!("accept failed: {e}")))?;
        served += 1;
        let reply = match handle_job(&mut stream, resolve.as_ref()) {
            Ok(body) => body,
            Err(e) => {
                let mut w = Writer::default();
                w.u8(REPLY_ERROR);
                w.str(&e.to_string());
                w.finish()
            }
        };
        // A client that disconnected mid-reply has nothing left to receive it.
        let _ = write_frame(&mut stream, &reply);
    }
}

/// Splits `rows` into `parts` contiguous ranges that differ in length by at most one.
fn split_rows(rows: usize, parts: usize) -> Vec<Range<usize>> {
    let base = rows / parts;
    let extra = rows % parts;
    let mut start = 0;
    (0..parts)
        .map(|i| {
            let len = base + usize::from(i < extra);
            let range = start..start + len;
            start += len;
            range
        })
        .collect()
}

enum Reply {
    Frame(DataFrame),
    Partial(Partial),
}

fn encode_job(source_name: &str, rows: Range<usize>, plan: &DistributedPlan) -> Vec<u8> {
    let mut w = Writer::default();
    w.u8(MSG_JOB);
    w.u8(PROTOCOL_VERSION);
    w.str(source_name);
    w.u64(rows.start as u64);
    w.u64(rows.end as u64);
    put_plan(&mut w, plan);
    w.finish()
}

fn put_plan(w: &mut Writer, plan: &DistributedPlan) {
    w.len(plan.stages.len());
    for stage in &plan.stages {
        match stage {
            Stage::Filter(e) => {
                w.u8(0);
                codec::put_expr(w, e);
            }
            Stage::Project(es) => {
                w.u8(1);
                w.len(es.len());
                for e in es {
                    codec::put_expr(w, e);
                }
            }
            Stage::WithColumn(e) => {
                w.u8(2);
                codec::put_expr(w, e);
            }
        }
    }
    match &plan.aggregate {
        None => w.u8(0),
        Some((keys, aggs)) => {
            w.u8(1);
            w.len(keys.len());
            for k in keys {
                w.str(k);
            }
            w.len(aggs.len());
            for a in aggs {
                codec::put_agg(w, a);
            }
        }
    }
}

fn get_plan(r: &mut Reader<'_>) -> Result<DistributedPlan, LazyError> {
    let nstages = r.count(1)?;
    let mut stages = Vec::with_capacity(nstages);
    for _ in 0..nstages {
        stages.push(match r.u8()? {
            0 => Stage::Filter(codec::get_expr(r)?),
            1 => {
                let n = r.count(2)?;
                Stage::Project((0..n).map(|_| codec::get_expr(r)).collect::<Result<_, _>>()?)
            }
            2 => Stage::WithColumn(codec::get_expr(r)?),
            other => return Err(LazyError::Invalid(format!("wire: unknown stage tag {other}"))),
        });
    }
    let aggregate = match r.u8()? {
        0 => None,
        1 => {
            let nkeys = r.count(4)?;
            let keys = (0..nkeys).map(|_| r.str()).collect::<Result<Vec<_>, _>>()?;
            let naggs = r.count(3)?;
            let aggs = (0..naggs).map(|_| codec::get_agg(r)).collect::<Result<Vec<_>, _>>()?;
            Some((keys, aggs))
        }
        other => return Err(LazyError::Invalid(format!("wire: bad aggregate flag {other}"))),
    };
    Ok(DistributedPlan { stages, aggregate })
}

/// Worker side: runs one job and returns the encoded reply.
fn handle_job(stream: &mut TcpStream, resolve: &ResolveFn) -> Result<Vec<u8>, LazyError> {
    let msg = read_frame(stream)?;
    let mut r = Reader::new(&msg);
    if r.u8()? != MSG_JOB {
        return Err(LazyError::Invalid("wire: expected a job message".into()));
    }
    if r.u8()? != PROTOCOL_VERSION {
        return Err(LazyError::Invalid("wire: protocol version mismatch".into()));
    }
    let source_name = r.str()?;
    let start =
        usize::try_from(r.u64()?).map_err(|_| LazyError::Invalid("wire: bad range".into()))?;
    let end =
        usize::try_from(r.u64()?).map_err(|_| LazyError::Invalid("wire: bad range".into()))?;
    if start > end {
        return Err(LazyError::Invalid("wire: range start after end".into()));
    }
    let plan = get_plan(&mut r)?;
    r.finish()?;

    let source = resolve(&source_name)?;
    let projection: Vec<String> = source.schema().into_iter().map(|f| f.name).collect();
    let rows = source.read(start..end, &projection)?;
    let local = plan.to_lazy(Arc::new(rows)).plan().clone();

    let mut w = Writer::default();
    if plan.aggregate.is_some() {
        let partial = exec::partial(&local, &CollectOptions::default())?;
        w.u8(REPLY_PARTIAL);
        codec::put_partial(&mut w, &partial);
    } else {
        let df = exec::execute(&local, &CollectOptions::default())?;
        w.u8(REPLY_FRAME);
        codec::put_frame(&mut w, &df)?;
    }
    Ok(w.finish())
}

/// Coordinator side: sends one job to `addr` and decodes the reply.
fn send_job(addr: SocketAddr, job: &[u8]) -> Result<Reply, LazyError> {
    let fail = |what: &str, e: &dyn std::fmt::Display| {
        LazyError::Invalid(format!("worker {addr}: {what}: {e}"))
    };
    let mut stream =
        TcpStream::connect_timeout(&addr, CONNECT_TIMEOUT).map_err(|e| fail("connect", &e))?;
    write_frame(&mut stream, job).map_err(|e| fail("send", &e))?;
    let reply = read_frame(&mut stream).map_err(|e| fail("receive", &e))?;
    let mut r = Reader::new(&reply);
    match r.u8()? {
        REPLY_FRAME => {
            let df = codec::get_frame(&mut r)?;
            r.finish()?;
            Ok(Reply::Frame(df))
        }
        REPLY_PARTIAL => {
            let p = codec::get_partial(&mut r)?;
            r.finish()?;
            Ok(Reply::Partial(p))
        }
        REPLY_ERROR => Err(LazyError::Invalid(format!("worker {addr}: {}", r.str()?))),
        other => Err(LazyError::Invalid(format!("worker {addr}: unknown reply tag {other}"))),
    }
}

fn write_frame(stream: &mut TcpStream, payload: &[u8]) -> std::io::Result<()> {
    let len = u32::try_from(payload.len())
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidInput, "message too large"))?;
    stream.write_all(&len.to_le_bytes())?;
    stream.write_all(payload)?;
    stream.flush()
}

fn read_frame(stream: &mut TcpStream) -> Result<Vec<u8>, LazyError> {
    let mut len = [0u8; 4];
    stream.read_exact(&mut len).map_err(|e| LazyError::Invalid(format!("wire: {e}")))?;
    let len = u32::from_le_bytes(len) as usize;
    if len > MAX_MESSAGE {
        return Err(LazyError::Invalid("wire: message exceeds size limit".into()));
    }
    let mut buf = vec![0u8; len];
    stream.read_exact(&mut buf).map_err(|e| LazyError::Invalid(format!("wire: {e}")))?;
    Ok(buf)
}
