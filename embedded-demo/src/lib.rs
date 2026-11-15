//! Fixed-point event-stream filter for a DAQ front end.
//!
//! This library models the part of a detector readout that runs on a
//! microcontroller or in FPGA-adjacent firmware: take raw ADC samples, convert
//! them to Q16.16 energies, and keep only events inside an energy window that
//! are not within a dead time of the previously accepted event. It uses no
//! floating point and no heap, so it is suitable for `no_std` targets without an
//! FPU.
//!
//! The crate depends on `tpt-dynamics-core` with default features off. The
//! library builds for `thumbv7em-none-eabihf` with `cargo build --lib` and does
//! not declare a panic handler (a firmware binary would supply one).
//!
//! # Example
//!
//! ```
//! use embedded_demo::{Filter, FilterConfig, q16_from_int};
//!
//! let cfg = FilterConfig {
//!     gain_q16: 1 << 16, // 1 ADC count = 1 unit
//!     lo_q16: q16_from_int(100), // accept energies in [100, 500)
//!     hi_q16: q16_from_int(500),
//!     deadtime_us: 10,
//! };
//! let mut filter = Filter::new(cfg);
//! let raw = [(0u32, 50u16), (5, 150), (30, 200), (31, 300)];
//! let kept = raw.iter().filter_map(|&(t, adc)| filter.push(t, adc)).count();
//! assert_eq!(kept, 2); // 50 is below the window; 31 is within dead time of 30
//! ```
#![no_std]
#![forbid(unsafe_code)]

/// Speed of light from `tpt-dynamics-core`, re-exported to show the core
/// crate is linked in its `no_std` configuration.
pub const SPEED_OF_LIGHT_M_S: f64 = tpt_dynamics_core::constants::SPEED_OF_LIGHT;

/// Q16.16 fixed-point value: the integer `v` represents `v / 65536`.
pub type Q16 = i32;

/// Converts an integer to Q16.16, saturating at the `i32` range.
#[must_use]
pub const fn q16_from_int(v: i32) -> Q16 {
    v.saturating_mul(1 << 16)
}

/// Multiplies two Q16.16 values with a 64-bit intermediate, saturating on overflow.
#[must_use]
pub fn q16_mul(a: Q16, b: Q16) -> Q16 {
    let wide = (i64::from(a) * i64::from(b)) >> 16;
    wide.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

/// One accepted event.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Event {
    /// Timestamp in microseconds since the start of the run.
    pub t_us: u32,
    /// Calibrated energy in Q16.16.
    pub energy_q16: Q16,
}

/// Static configuration of the filter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FilterConfig {
    /// Gain applied to raw ADC counts, in Q16.16.
    pub gain_q16: Q16,
    /// Inclusive lower energy bound, Q16.16.
    pub lo_q16: Q16,
    /// Exclusive upper energy bound, Q16.16.
    pub hi_q16: Q16,
    /// Minimum separation between accepted events, in microseconds.
    pub deadtime_us: u32,
}

/// Running counters kept by the filter.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counters {
    /// Samples presented to the filter.
    pub seen: u32,
    /// Samples rejected by the energy window.
    pub outside_window: u32,
    /// Samples rejected by the dead time.
    pub dead: u32,
    /// Samples accepted.
    pub accepted: u32,
}

/// Streaming energy-window and dead-time filter.
#[derive(Clone, Copy, Debug)]
pub struct Filter {
    cfg: FilterConfig,
    last_accept_us: Option<u32>,
    counters: Counters,
}

impl Filter {
    /// Create a filter with the given configuration.
    #[must_use]
    pub fn new(cfg: FilterConfig) -> Self {
        Self { cfg, last_accept_us: None, counters: Counters::default() }
    }

    /// Current counters.
    #[must_use]
    pub fn counters(&self) -> Counters {
        self.counters
    }

    /// Convert a raw ADC count to an energy and test it. Returns the event if it
    /// is accepted.
    pub fn push(&mut self, t_us: u32, adc: u16) -> Option<Event> {
        let energy = q16_mul(q16_from_int(i32::from(adc)), self.cfg.gain_q16);
        self.push_energy(t_us, energy)
    }

    /// Test an event whose energy is already in Q16.16.
    pub fn push_energy(&mut self, t_us: u32, energy_q16: Q16) -> Option<Event> {
        self.counters.seen = self.counters.seen.saturating_add(1);
        if energy_q16 < self.cfg.lo_q16 || energy_q16 >= self.cfg.hi_q16 {
            self.counters.outside_window = self.counters.outside_window.saturating_add(1);
            return None;
        }
        if let Some(last) = self.last_accept_us {
            // Wrapping subtraction handles the u32 timestamp rollover.
            if t_us.wrapping_sub(last) < self.cfg.deadtime_us {
                self.counters.dead = self.counters.dead.saturating_add(1);
                return None;
            }
        }
        self.last_accept_us = Some(t_us);
        self.counters.accepted = self.counters.accepted.saturating_add(1);
        Some(Event { t_us, energy_q16 })
    }

    /// Filter a block of raw samples into a fixed-size output buffer.
    ///
    /// Returns the number of events written. Accepted events beyond the length of
    /// `out` are counted in [`Counters::accepted`] but not stored.
    pub fn process_block(&mut self, samples: &[(u32, u16)], out: &mut [Event]) -> usize {
        let mut n = 0;
        for &(t, adc) in samples {
            if let Some(ev) = self.push(t, adc) {
                if n < out.len() {
                    out[n] = ev;
                    n += 1;
                }
            }
        }
        n
    }
}
