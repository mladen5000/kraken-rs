use std::collections::HashSet;
use std::iter::FromIterator;
use std::ops::AddAssign;

pub struct ReadCounts<K> {
    n_reads: u64,
    n_kmers: u64,
    kmers: K,
}

impl<K> ReadCounts<K>
where
    K: FromIterator<u64> + IntoIterator<Item = u64>,
{
    pub fn new() -> Self {
        Self {
            n_reads: 0,
            n_kmers: 0,
            kmers: K::from_iter(std::iter::empty()),
        }
    }

    pub fn with_values(n_reads: u64, n_kmers: u64) -> Self {
        Self {
            n_reads,
            n_kmers,
            kmers: K::from_iter(std::iter::empty()),
        }
    }

    pub fn read_count(&self) -> u64 {
        self.n_reads
    }

    pub fn increment_read_count(&mut self) {
        self.n_reads += 1;
    }

    pub fn kmer_count(&self) -> u64 {
        self.n_kmers
    }

    pub fn distinct_kmer_count(&self) -> usize {
        self.kmers.into_iter().collect::<HashSet<_>>().len()
    }

    pub fn add_kmer(&mut self, kmer: u64) {
        self.n_kmers += 1;
        self.kmers = K::from_iter(self.kmers.into_iter().chain(std::iter::once(kmer)));
    }
}

impl<K> AddAssign for ReadCounts<K>
where
    K: FromIterator<u64> + IntoIterator<Item = u64>,
{
    fn add_assign(&mut self, other: Self) {
        self.n_reads += other.n_reads;
        self.n_kmers += other.n_kmers;
        self.kmers = K::from_iter(self.kmers.into_iter().chain(other.kmers.into_iter()));
    }
}

impl<K> PartialEq for ReadCounts<K>
where
    K: FromIterator<u64> + IntoIterator<Item = u64>,
{
    fn eq(&self, other: &Self) -> bool {
        self.n_reads == other.n_reads && self.n_kmers == other.n_kmers
    }
}

impl<K> PartialOrd for ReadCounts<K>
where
    K: FromIterator<u64> + IntoIterator<Item = u64>,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.n_reads < other.n_reads {
            Some(std::cmp::Ordering::Less)
        } else if self.n_reads == other.n_reads && self.n_kmers < other.n_kmers {
            Some(std::cmp::Ordering::Less)
        } else {
            Some(std::cmp::Ordering::Greater)
        }
    }
}
