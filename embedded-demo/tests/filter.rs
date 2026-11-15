//! Host-run tests for the fixed-point filter. The test harness uses std; the
//! library itself is no_std.

use embedded_demo::{q16_from_int, q16_mul, Event, Filter, FilterConfig, SPEED_OF_LIGHT_M_S};

fn cfg() -> FilterConfig {
    FilterConfig {
        gain_q16: 1 << 16,
        lo_q16: q16_from_int(100),
        hi_q16: q16_from_int(500),
        deadtime_us: 10,
    }
}

#[test]
fn window_and_deadtime() {
    let mut f = Filter::new(cfg());
    assert_eq!(f.push(0, 99), None);
    assert_eq!(f.push(1, 100), Some(Event { t_us: 1, energy_q16: q16_from_int(100) }));
    assert_eq!(f.push(5, 200), None); // inside dead time
    assert!(f.push(11, 200).is_some()); // exactly deadtime after the last accept
    assert_eq!(f.push(30, 500), None); // upper bound is exclusive
    let c = f.counters();
    assert_eq!((c.seen, c.accepted, c.dead, c.outside_window), (5, 2, 1, 2));
}

#[test]
fn timestamp_rollover_is_handled() {
    let mut f = Filter::new(cfg());
    assert!(f.push(u32::MAX - 2, 200).is_some());
    assert_eq!(f.push(u32::MAX, 200), None); // 2 us later: dead
    assert!(f.push(7, 200).is_some()); // 10 us later across the wrap
}

#[test]
fn q16_arithmetic() {
    assert_eq!(q16_mul(q16_from_int(3), q16_from_int(4)), q16_from_int(12));
    assert_eq!(q16_mul(i32::MAX, i32::MAX), i32::MAX); // saturates
    assert_eq!(q16_from_int(i32::MAX), i32::MAX);
    assert_eq!(q16_from_int(-1), -65536);
}

#[test]
fn gain_scales_adc_counts() {
    let mut c = cfg();
    c.gain_q16 = 1 << 15; // 0.5 energy units per ADC count
    let mut f = Filter::new(c);
    assert_eq!(f.push(0, 199), None); // 99.5 < 100
    assert!(f.push(100, 200).is_some()); // 100.0 is accepted
}

#[test]
fn block_processing_matches_streaming() {
    let samples: Vec<(u32, u16)> = (0..500u32).map(|i| (i * 3, ((i * 37) % 700) as u16)).collect();
    let mut streaming = Filter::new(cfg());
    let expected: Vec<Event> = samples.iter().filter_map(|&(t, a)| streaming.push(t, a)).collect();

    let mut block = Filter::new(cfg());
    let mut out = [Event { t_us: 0, energy_q16: 0 }; 512];
    let n = block.process_block(&samples, &mut out);
    assert_eq!(&out[..n], &expected[..]);
    assert_eq!(block.counters(), streaming.counters());
}

#[test]
fn block_output_is_bounded() {
    let samples: Vec<(u32, u16)> = (0..100u32).map(|i| (i * 100, 200)).collect();
    let mut f = Filter::new(cfg());
    let mut out = [Event { t_us: 0, energy_q16: 0 }; 10];
    assert_eq!(f.process_block(&samples, &mut out), 10);
    assert_eq!(f.counters().accepted, 100);
}

#[test]
fn core_constant_is_linked() {
    assert!((SPEED_OF_LIGHT_M_S - 299_792_458.0).abs() < 1.0);
}
