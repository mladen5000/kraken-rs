use anyhow::Error;

pub enum ErrorKind {
    IoError,
    FormatError,
    ParseError,
    DatabaseError,
    TaxonomyError,
}

pub fn exit_with_error(kind: ErrorKind, message: &str) -> ! {
    eprintln!("Error: {}", message);
    std::process::exit(1);
}
