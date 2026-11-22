# Changelog

All notable changes to `embedded-demo` are documented in this file. The workspace
changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md). This crate is a
standalone demonstration and is not published (`publish = false`).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `Filter`: a streaming energy-window and dead-time filter over raw ADC samples, in
  Q16.16 fixed point. Events are accepted inside `[lo, hi)` and outside the dead time
  of the previous accepted event.
- `Counters`: `seen`, `outside_window`, `dead` and `accepted`, all saturating.
- `Filter::process_block`: filters a block of samples into a caller-supplied buffer.
  Accepted events beyond the buffer are counted, but not stored.
- `q16_from_int` and `q16_mul` fixed-point helpers.
- Host tests in `tests/filter.rs`: window and dead-time rules, timestamp rollover,
  Q16.16 arithmetic, gain scaling, block and streaming equivalence, bounded block
  output, and linkage of the `tpt-dynamics-core` constant.
- Continuous integration builds the crate for `thumbv7em-none-eabihf` on every change.
