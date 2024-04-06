/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file was originally developed for the KrakenUniq taxonomic classification system.
 */

// Note: This class provides counter information for read and k-mer data
// as well as allowing counting/estimation of distinct k-mers via a container
// type that is passed in.

use std::collections::HashSet;

use crate::hyperloglogplus::{murmurhash3_finalizer, HyperLogLogPlusMinus};

pub struct ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    n_reads: u64,
    n_kmers: u64,
    kmers: CONTAINER,
}

impl<CONTAINER> ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    pub fn distinct_kmer_count(&self) -> u64 {
        self.kmers.size()
    }

    pub fn new(kmers: CONTAINER) -> Self {
        Self {
            n_reads: 0,
            n_kmers: 0,
            kmers,
        }
    }

    pub fn with_counts(n_reads: u64, n_kmers: u64, kmers: CONTAINER) -> Self {
        Self {
            n_reads,
            n_kmers,
            kmers,
        }
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers.insert(kmer);
    }
    pub fn add(&mut self, other: &ReadCounts<CONTAINER>) {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        for kmer in other.kmers.iter() {
            self.kmers.insert(kmer);
        }
    }
}

impl<CONTAINER> PartialEq for ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    fn eq(&self, other: &Self) -> bool {
        self.n_reads == other.n_reads && self.n_kmers == other.n_kmers
    }
}

impl<CONTAINER> Eq for ReadCounts<CONTAINER> where CONTAINER: KmerContainer {}

impl<CONTAINER> PartialOrd for ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<CONTAINER> Ord for ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.n_reads.cmp(&other.n_reads) {
            std::cmp::Ordering::Equal => self.n_kmers.cmp(&other.n_kmers),
            ord => ord,
        }
    }
}

pub trait KmerContainer {
    fn new() -> Self;
    fn insert(&mut self, kmer: u64);
    fn size(&self) -> u64;
    fn iter(&self) -> Box<dyn Iterator<Item = u64> + '_>;
}

impl KmerContainer for HashSet<u64> {
    fn new() -> Self {
        HashSet::new()
    }

    fn insert(&mut self, kmer: u64) {
        self.insert(kmer);
    }

    fn size(&self) -> u64 {
        self.len() as u64
    }
    fn iter(&self) -> Box<dyn Iterator<Item = u64> + '_> {
        Box::new(self.iter().copied())
    }
}

impl KmerContainer for HyperLogLogPlusMinus {
    fn new() -> Self {
        HyperLogLogPlusMinus::new(12, true, murmurhash3_finalizer)
    }

    fn insert(&mut self, kmer: u64) {
        self.insert(kmer);
    }

    fn size(&self) -> u64 {
        self.cardinality()
    }
    fn iter(&self) -> Box<dyn Iterator<Item = u64> + '_> {
        Box::new(self.sparse_list.iter().map(|&kmer| kmer as u64))
    }
}

impl<CONTAINER> Default for ReadCounts<CONTAINER>
where
    CONTAINER: KmerContainer,
{
    fn default() -> Self {
        Self {
            n_reads: 0,
            n_kmers: 0,
            kmers: CONTAINER::new(),
        }
    }
}
