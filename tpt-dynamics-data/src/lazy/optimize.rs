//! Rule-based plan optimiser: constant folding, predicate pushdown and projection pushdown.
//!
//! Every rule preserves results exactly. Predicate pushdown only moves a filter
//! past operators whose outputs it does not depend on. Projection pushdown
//! removes unused columns and unused with-column nodes, and tells the scan to
//! read only the columns the plan needs.

use std::collections::BTreeSet;

use super::expr::{fold_constants, Expr};
use super::plan::Plan;

/// Applies all rules: push filters down first, then prune columns.
pub(crate) fn optimize(plan: Plan) -> Plan {
    let plan = fold_plan(plan);
    let plan = pushdown_predicates(plan);
    prune(plan, None)
}

/// Folds literal-only operators inside every expression of the plan. Aggregate
/// inputs are left as written.
fn fold_plan(plan: Plan) -> Plan {
    match plan {
        Plan::Scan { .. } => plan,
        Plan::Filter { input, predicate } => Plan::Filter {
            input: Box::new(fold_plan(*input)),
            predicate: fold_constants(predicate),
        },
        Plan::Project { input, exprs } => Plan::Project {
            input: Box::new(fold_plan(*input)),
            exprs: exprs.into_iter().map(fold_constants).collect(),
        },
        Plan::WithColumn { input, expr } => {
            Plan::WithColumn { input: Box::new(fold_plan(*input)), expr: fold_constants(expr) }
        }
        Plan::Aggregate { input, keys, aggs } => {
            Plan::Aggregate { input: Box::new(fold_plan(*input)), keys, aggs }
        }
    }
}

fn pushdown_predicates(plan: Plan) -> Plan {
    match plan {
        Plan::Scan { .. } => plan,
        Plan::Filter { input, predicate } => push_filter(pushdown_predicates(*input), predicate),
        Plan::Project { input, exprs } => {
            Plan::Project { input: Box::new(pushdown_predicates(*input)), exprs }
        }
        Plan::WithColumn { input, expr } => {
            Plan::WithColumn { input: Box::new(pushdown_predicates(*input)), expr }
        }
        Plan::Aggregate { input, keys, aggs } => {
            Plan::Aggregate { input: Box::new(pushdown_predicates(*input)), keys, aggs }
        }
    }
}

/// `true` when `expr` outputs column `name` unchanged (`col(name)` or
/// `col(name).alias(name)`).
fn is_identity(expr: &Expr, name: &str) -> bool {
    expr.output_name() == Some(name) && matches!(expr.unaliased(), Expr::Col(n) if n == name)
}

/// Places `predicate` as low in the plan as its column references allow.
fn push_filter(input: Plan, predicate: Expr) -> Plan {
    let refs = predicate.referenced_columns();
    match input {
        // Only identity outputs can be read from below the projection.
        Plan::Project { input: inner, exprs }
            if refs.iter().all(|c| exprs.iter().any(|e| is_identity(e, c))) =>
        {
            Plan::Project { input: Box::new(push_filter(*inner, predicate)), exprs }
        }
        // The predicate must not read the column the with-column overwrites.
        Plan::WithColumn { input: inner, expr }
            if expr.output_name().is_some_and(|n| !refs.contains(n)) =>
        {
            Plan::WithColumn { input: Box::new(push_filter(*inner, predicate)), expr }
        }
        // Filtering on group keys commutes with grouping.
        Plan::Aggregate { input: inner, keys, aggs } if refs.iter().all(|c| keys.contains(c)) => {
            Plan::Aggregate { input: Box::new(push_filter(*inner, predicate)), keys, aggs }
        }
        other => Plan::Filter { input: Box::new(other), predicate },
    }
}

/// Removes unused outputs and tells scans which columns to read. `required` is
/// the set of column names the parent needs, or `None` for all of them.
fn prune(plan: Plan, required: Option<BTreeSet<String>>) -> Plan {
    match plan {
        Plan::Scan { source, projection } => {
            let projection = match &required {
                Some(req) => Some(
                    source
                        .schema()
                        .into_iter()
                        .map(|f| f.name)
                        .filter(|n| req.contains(n))
                        .collect(),
                ),
                None => projection,
            };
            Plan::Scan { source, projection }
        }
        Plan::Filter { input, predicate } => {
            let req = required.map(|mut r| {
                r.extend(predicate.referenced_columns());
                r
            });
            Plan::Filter { input: Box::new(prune(*input, req)), predicate }
        }
        Plan::Project { input, exprs } => {
            let kept: Vec<Expr> = exprs
                .into_iter()
                .filter(|e| match (&required, e.output_name()) {
                    (Some(req), Some(n)) => req.contains(n),
                    _ => true,
                })
                .collect();
            let refs: BTreeSet<String> = kept.iter().flat_map(Expr::referenced_columns).collect();
            Plan::Project { input: Box::new(prune(*input, Some(refs))), exprs: kept }
        }
        Plan::WithColumn { input, expr } => {
            let Some(name) = expr.output_name().map(str::to_string) else {
                return Plan::WithColumn { input, expr };
            };
            if required.as_ref().is_some_and(|req| !req.contains(&name)) {
                // Nobody reads the new column: drop the node entirely.
                return prune(*input, required);
            }
            let req = required.map(|mut r| {
                r.remove(&name);
                r.extend(expr.referenced_columns());
                r
            });
            Plan::WithColumn { input: Box::new(prune(*input, req)), expr }
        }
        Plan::Aggregate { input, keys, aggs } => {
            let mut req: BTreeSet<String> = keys.iter().cloned().collect();
            for a in &aggs {
                if let Some(e) = a.input() {
                    req.extend(e.referenced_columns());
                }
            }
            Plan::Aggregate { input: Box::new(prune(*input, Some(req))), keys, aggs }
        }
    }
}
