# KRAKEN-RS Development Guide

## Build Commands
- Build all: `./build.sh build`
- Install: `./build.sh install`
- Clean: `./build.sh clean`
- Build specific binary: `cargo build --bin <binary_name> --release`
- Run tests: `cargo test`
- Run single test: `cargo test <test_name>`
- Run specific module tests: `cargo test --package kraken-rs --lib <module_name>::tests`
- Lint with clippy: `cargo clippy`
- Build with features: `cargo build --features <feature_name>`

## Code Style Guidelines
- **Error Handling**: Use anyhow for error propagation; custom ErrorKind enum for specific errors; add context with `.context()`
- **Naming**: Use snake_case for variables/functions, CamelCase for types/structs, test_* for test functions
- **Module Structure**: Each functional component has its own module; tests in module-level `tests` submodules with `#[cfg(test)]`
- **Imports**: Group standard library, external crates, then internal modules
- **Documentation**: Use doc comments `///` for public API and functions; explain unsafe blocks
- **Threading**: Prefer rayon for parallel processing; use threadpool for custom thread management; ensure thread safety
- **Memory Management**: Use mmap for large files; implement Drop for cleanup; use unsafe with careful validation
- **Performance**: Optimize critical paths; use bit-packing and atomic operations where appropriate