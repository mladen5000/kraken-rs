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

## Code Style Guidelines
- **Error Handling**: Use anyhow for error propagation; custom ErrorKind enum for specific errors
- **Naming**: Use snake_case for variables/functions, CamelCase for types/structs
- **Module Structure**: Each functional component has its own module
- **Testing**: Unit tests in module-level `tests` submodules with `#[cfg(test)]`
- **Imports**: Group standard library, external crates, then internal modules
- **Documentation**: Use doc comments `///` for public API and functions
- **Threading**: Prefer rayon for parallel processing; use threadpool for custom thread management
- **Memory Management**: Utilize mmap for large file handling