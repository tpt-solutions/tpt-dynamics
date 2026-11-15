//! A minimal stderr logger for the `log` facade.

use log::{LevelFilter, Log, Metadata, Record};

struct StderrLogger;

static LOGGER: StderrLogger = StderrLogger;

impl Log for StderrLogger {
    fn enabled(&self, metadata: &Metadata<'_>) -> bool {
        metadata.level() <= log::max_level()
    }

    fn log(&self, record: &Record<'_>) {
        if self.enabled(record.metadata()) {
            eprintln!("[{}] {}", record.level(), record.args());
        }
    }

    fn flush(&self) {}
}

/// Installs the stderr logger at `level`. Calling it again is harmless.
pub(crate) fn init(level: LevelFilter) {
    let _ = log::set_logger(&LOGGER);
    log::set_max_level(level);
}
