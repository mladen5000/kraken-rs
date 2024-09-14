use clap::{App, Arg};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

fn main() {
    let matches = App::new("File Reader")
        .version("1.0")
        .author("Your Name")
        .about("Reads a file and prints its contents")
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("Sets the input file to read")
                .required(true)
                .takes_value(true),
        )
        .get_matches();

    let file_path = matches.value_of("input").unwrap();
    if let Err(err) = read_file(file_path) {
        eprintln!("Hello, world!");
    }
}

fn read_file(file_path: &str) -> io::Result<()> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        println!("{}", line?);
    }

    Ok(())
}
