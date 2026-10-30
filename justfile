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

# Everything CI runs
check: fmt-check clippy test deny

# Build tpt-chem-core for a bare-metal no_std target
no-std:
    cargo build -p tpt-chem-core --target thumbv6m-none-eabi --no-default-features

# Release build
release:
    cargo build --release

# Build rustdoc for the workspace, denying warnings
doc:
    RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --all-features

# Run every umbrella-crate example (each is feature-gated via required-features)
examples:
    cargo run -p tpt-chemistry --example h2_hf --features core,quantum --release
    cargo run -p tpt-chemistry --example water_hf --features core,quantum --release
    cargo run -p tpt-chemistry --example lj_argon_nve --features core,md --release
    cargo run -p tpt-chemistry --example nacl_pme --features core,md --release
    cargo run -p tpt-chemistry --example lotka_volterra_ssa --features core,kinetics --release
    cargo run -p tpt-chemistry --example michaelis_menten --features core,kinetics --release
    cargo run -p tpt-chemistry --example water_pme --features core,md --release
    cargo run -p tpt-chemistry --example silicon_xrd --features core,crystal --release
    cargo run -p tpt-chemistry --example nacl_xrd --features core,crystal --release
    cargo run -p tpt-chemistry --example cif_to_xrd --features core,io,crystal --release
    cargo run -p tpt-chemistry --example cpp_calcium_phosphate --features core,aqueous --release

# Run all benchmarks (none yet; the criterion job is still open in todo.md)
bench:
    cargo bench --workspace
