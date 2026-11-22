# embedded-demo

A fixed-point event-stream filter for a detector data-acquisition (DAQ) front end.
It shows how the pieces of [tpt-dynamics](../README.md) fit on a microcontroller or
in firmware next to an FPGA: raw ADC samples go in, and only the events that matter
come out.

This is a demonstration, not a production firmware image. It is included to prove
that the core crate builds for a bare-metal target and to show a realistic
`no_std` integration.

## What it does

For each raw sample `(time in µs, ADC count)` the filter:

1. Converts the ADC count to an energy in Q16.16 fixed point, using a gain.
2. Rejects the event if the energy is outside the window `[lo, hi)`.
3. Rejects the event if it falls within `deadtime_us` of the previous **accepted**
   event. Rejected events do not restart the dead time.
4. Otherwise accepts it, and records it as an `Event { t_us, energy_q16 }`.

Every decision is counted in `Counters`: `seen`, `outside_window`, `dead` and
`accepted`. All counters saturate rather than wrap.

Timestamps are `u32` microseconds and may roll over. The dead-time test uses a
wrapping subtraction, so a rollover does not block events.

## Properties

- **No floating point and no heap.** All arithmetic is `i32` Q16.16 (`Q16`).
  Multiplication goes through `q16_mul`. The filter suits a core without an FPU.
- **`no_std`, no `unsafe`.** The library is `#![no_std]` and `#![forbid(unsafe_code)]`.
- **Bounded output.** `process_block` writes accepted events into a caller-supplied
  slice. Events that do not fit are still counted in `Counters::accepted`, but not
  stored.
- **Depends on `tpt-dynamics-core` with default features off** (so not its `std`
  feature). The crate's `SPEED_OF_LIGHT_M_S` constant is taken from
  `tpt-dynamics-core::constants`.

## Example

```rust
use embedded_demo::{Filter, FilterConfig, q16_from_int};

let cfg = FilterConfig {
    gain_q16: 1 << 16,          // 1 ADC count = 1 unit
    lo_q16: q16_from_int(100),  // accept energies in [100, 500)
    hi_q16: q16_from_int(500),
    deadtime_us: 10,
};
let mut filter = Filter::new(cfg);
let raw = [(0u32, 50u16), (5, 150), (30, 200), (31, 300)];
let kept = raw.iter().filter_map(|&(t, adc)| filter.push(t, adc)).count();
assert_eq!(kept, 2); // 50 is below the window; 31 is within dead time of 30
```

## Building

This crate is a **separate workspace** with its own `Cargo.lock`. It is not a member of
the root workspace, so `cargo build` at the repository root does not build it. Run
commands from this directory:

```text
rustup target add thumbv7em-none-eabihf
cargo build --locked --target thumbv7em-none-eabihf
```

The library builds for the target without a panic handler. A firmware binary would
supply one. CI runs the same command on every change (the `embedded-demo` job in
`.github/workflows/ci.yml`).

Host tests run with:

```text
cargo test
```

The tests in `tests/filter.rs` cover the window and dead-time rules, timestamp
rollover, Q16.16 arithmetic, gain scaling, block and streaming equivalence, bounded
block output, and linkage of the core constant.

## Not included

- Hardware drivers, interrupts, a clock source or a panic handler.
- Any interface to a real ADC or FPGA. Samples are supplied by the caller.
- Energy calibration. The gain is a single fixed-point factor.

## Licence

Licensed under either of [MIT](../LICENSE-MIT) or [Apache-2.0](../LICENSE-APACHE), at
your option. Part of [tpt-dynamics](../README.md). Not published to crates.io.
