use std::collections::{HashMap, HashSet};
use std::error::Error;

use flate2::read::GzDecoder;
use futures::stream::{FuturesUnordered, StreamExt};
use log;
use reqwest;
use tokio;

#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        println!("Usage: {} <input_file>", args[0]);
        return;
    }
    let input_file = &args[1];

    let file = tokio::fs::File::open(input_file).unwrap();
    let reader = GzDecoder::new(file);
    let mut reader = tokio::io::BufReader::new(reader);
    let mut reader = futures::io::AsyncReadExt::reader(reader);
    let mut buf = Vec::new();
    let mut all = HashSet::new();
    let mut futures = FuturesUnordered::new();

    loop {
        match reader.read_until(b'\n', &mut buf).await {
            Ok(n) => {
                if n == 0 {
                    break;
                }
                let line = String::from_utf8_lossy(&buf[..n]);
                let mut parts = line.split('\t');
                let taxid = parts.next().unwrap();
                let name = parts.next().unwrap();
                let taxid = taxid.parse::<u32>().unwrap();
                let name = name.to_string();
                futures.push(tokio::spawn(async move {
                    let res = reqwest::get(&format!(
                        "http://api.ncbi.nlm.nih.gov/taxonomy/lookup?db=nucl&id={}",
                        taxid
                    ))
                    .await
                    .unwrap();
                    let body = res.text().await.unwrap();
                    let mut parts = body.split('\t');
                    let taxid = parts.next().unwrap();
                    let name = parts.next().unwrap();
                    let taxid = taxid.parse::<u32>().unwrap();
                    let name = name.to_string();
                    (taxid, name)
                }));
                buf.clear();
            }
            Err(e) => {
                log::error!("Error reading file: {}", e);
                break;
            }
        }
    }
    while let Some(res) = futures.next().await {
        let (taxid, name) = res.unwrap();
        all.insert((taxid, name));
    }
    for (taxid, name) in &all {
        println!("{} {}", taxid, name);
    }
}
