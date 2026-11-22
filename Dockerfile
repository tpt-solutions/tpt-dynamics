# Multi-stage build of the tptdyn command-line tool.
#
#   docker build -t tptdyn .
#   docker run --rm -v "$PWD:/work" -w /work tptdyn --help
#
# The builder uses the MSRV from Cargo.toml (workspace.package.rust-version). The
# runtime image holds only the binary and the licence files.

ARG RUST_VERSION=1.95
ARG DEBIAN_VERSION=bookworm

FROM rust:${RUST_VERSION}-${DEBIAN_VERSION} AS builder
WORKDIR /src
COPY . .
# --locked uses the committed Cargo.lock. vendor/ is part of the build (the
# [patch.crates-io] entry for target-lexicon points at it), so it is not excluded.
RUN cargo build --release --locked -p tpt-dynamics-cli \
    && cp target/release/tptdyn /usr/local/bin/tptdyn

FROM debian:bookworm-slim AS runtime
RUN useradd --system --create-home --uid 10001 tptdyn
COPY --from=builder /usr/local/bin/tptdyn /usr/local/bin/tptdyn
COPY LICENSE-MIT LICENSE-APACHE NOTICE /usr/share/doc/tptdyn/
USER tptdyn
WORKDIR /work
ENTRYPOINT ["tptdyn"]
CMD ["--help"]
