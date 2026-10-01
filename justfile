default:
  @just --list

# Format all code
fmt:
    cargo fmt --all

# Check formatting without writing
fmt-check:
    cargo fmt --all -- --check

# Lint everything with clippy, denying warnings
clippy:
    cargo clippy --all-targets --all-features -- -D warnings

# Run the full test suite
test:
    cargo test --workspace

# License / ban / source audit (no Apache-2.0-only deps allowed)
deny:
    cargo deny check

# Everything CI runs, minus Kani
check: fmt-check clippy test deny

# Build tpt-chem-core for a bare-metal no_std target
no-std:
    cargo build -p tpt-chem-core --target thumbv6m-none-eabi --no-default-features

# Bounded model checking (Kani does not support Windows; run on Linux/WSL)
kani:
    cargo kani --workspace

# Release build
release:
    cargo build --release
