use crate::kv_store::murmurhash3;

use getopts::Options;
use memmap2::MmapMut;
use rayon::prelude::*; // Assuming memmap2 crate is used for memory mapping

use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::{Arc, Mutex};

const LOCK_ZONES: usize = 256;

pub struct CompactHashCell {
    data: u32,
}

impl CompactHashCell {
    fn hashed_key(&self, value_bits: usize) -> u64 {
        (self.data as u64) >> value_bits
    }

    fn value(&self, value_bits: usize) -> u64 {
        self.data as u64 & ((1 << value_bits) - 1)
    }

    fn populate(
        &mut self,
        compacted_key: u64,
        val: u64,
        key_bits: usize,
        value_bits: usize,
    ) -> Result<(), Box<dyn Error>> {
        if key_bits + value_bits != 32 {
            return Err("key len and value len don't sum to 32".into());
        }
        if key_bits == 0 || value_bits == 0 {
            return Err("key len and value len must be nonzero".into());
        }
        let max_value = (1 << value_bits) - 1;
        if max_value < val {
            return Err("value len too small for value".into());
        }
        self.data = (compacted_key << value_bits) as u32;
        self.data |= val as u32;
        Ok(())
    }
}

struct CompactHashTable {
    capacity: usize,
    size: AtomicUsize,
    key_bits: usize,
    value_bits: usize,
    table: Vec<Arc<Mutex<CompactHashCell>>>, // Use Arc for shared ownership across threads
    file_backed: bool,
    locks: Vec<Mutex<()>>, // Replacing omp_lock_t with Mutex
    mmap: Option<MmapMut>, // Optional memory-mapped file backing
}

impl CompactHashTable {
    pub fn new(capacity: usize, key_bits: usize, value_bits: usize) -> io::Result<Self> {
        let table = (0..capacity)
            .map(|_| Arc::new(Mutex::new(CompactHashCell { data: 0 })))
            .collect();
        let locks = (0..LOCK_ZONES).map(|_| Mutex::new(())).collect();

        Ok(Self {
            capacity,
            size: AtomicUsize::new(0),
            key_bits,
            value_bits,
            table,
            file_backed: false,
            locks,
            mmap: None,
        })
    }

    pub fn with_file<P: AsRef<Path>>(
        path: P,
        capacity: usize,
        key_bits: usize,
        value_bits: usize,
    ) -> io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapMut::map_mut(&file)? };
        // Further setup based on the file's contents

        Ok(Self {
            capacity,
            size: AtomicUsize::new(0), // Change the type to AtomicUsize
            key_bits,
            value_bits,
            table: vec![], // Initialize properly based on mmap and capacity
            file_backed: true,
            locks: (0..LOCK_ZONES).map(|_| Mutex::new(())).collect(),
            mmap: Some(mmap),
        })
    }

    pub fn get(&self, key: u64) -> Option<u64> {
        let hc = murmurhash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step = 0;

        loop {
            let cell = self.table[idx].lock().unwrap();
            if cell.value(self.value_bits) == 0 {
                // Value of 0 means data is 0, saves work
                break None;
            }
            if cell.hashed_key(self.value_bits) == compacted_key {
                return Some(cell.value(self.value_bits));
            }
            if step == 0 {
                step = self.second_hash(hc) as usize;
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break None; // Search over, we've exhausted the table
            }
        }
    }
    pub fn find_index(&self, key: u64) -> Option<usize> {
        let hc = murmurhash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = (hc % self.capacity as u64) as usize;
        let first_idx = idx;
        let mut step = 0;

        loop {
            let cell = self.table[idx].lock().unwrap();
            if cell.value(self.value_bits) == 0 {
                return None;
            }
            if cell.hashed_key(self.value_bits) == compacted_key {
                return Some(idx);
            }
            if step == 0 {
                step = self.second_hash(hc) as usize;
            }
            idx = (idx + step) % self.capacity;
            if idx == first_idx {
                break None;
            }
        }
    }
    pub fn compare_and_set(&self, key: u64, new_value: u64, old_value: &mut u64) -> bool {
        if self.file_backed {
            return false;
        }
        if new_value == 0 {
            return false;
        }

        let hc = murmurhash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = (hc % self.capacity as u64) as usize;
        let mut set_successful = false;

        let zone = idx % LOCK_ZONES;
        let _lock = self.locks[zone].lock().unwrap(); // Lock the corresponding zone for thread safety

        {
            let cell = &mut *self.table[idx].lock().unwrap(); // Lock the cell for exclusive access

            if *old_value == cell.value(self.value_bits) {
                cell.populate(compacted_key, new_value, self.key_bits, self.value_bits);
                if *old_value == 0 {
                    self.size.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                    // Safely increment size, assuming `size` is an `AtomicUsize`
                }
                *old_value = cell.value(self.value_bits); // Update the old value to the new value
                set_successful = true;
            } else {
                *old_value = cell.value(self.value_bits); // Update old_value if the comparison fails
            }
        }

        set_successful
    }

    pub fn direct_compare_and_set(
        &self,
        idx: usize,
        key: u64,
        new_value: u64,
        old_value: &mut u64,
    ) -> bool {
        let hc = murmurhash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let zone = idx % LOCK_ZONES;
        let _lock = self.locks[zone].lock().unwrap(); // Lock the corresponding zone

        let mut cell = self.table[idx].lock().unwrap();
        if *old_value == cell.value(self.value_bits) {
            cell.populate(compacted_key, new_value, self.key_bits, self.value_bits);
            if *old_value == 0 {
                // Assuming size is an AtomicUsize
                self.size.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            }
            *old_value = cell.value(self.value_bits);
            true
        } else {
            *old_value = cell.value(self.value_bits);
            false
        }
    }

    pub fn get_value_counts(&self) -> HashMap<u64, usize> {
        let value_counts: HashMap<u64, usize> = self
            .table
            .par_iter()
            .map(|cell| {
                let cell = cell.lock().unwrap();
                let value = cell.value(self.value_bits);
                if value != 0 {
                    Some((value, 1))
                } else {
                    None
                }
            })
            .filter_map(|x| x)
            .collect::<Vec<_>>()
            .into_iter()
            .fold(HashMap::new(), |mut acc, (val, count)| {
                *acc.entry(val).or_insert(0) += count;
                acc
            });
        value_counts
    }
    fn second_hash(&self, first_hash: u64) -> u64 {
        #[cfg(feature = "linear_probing")]
        let hash_value = 1;

        #[cfg(not(feature = "linear_probing"))]
        let hash_value = (first_hash >> 8) | 1;

        hash_value
    }
}

struct CmdOptions {
    hashtable_filename: String,
    taxonomy_filename: String,
    options_filename: String,
    output_filename: String,
    use_mpa_style: bool,
    report_zeros: bool,
    skip_counts: bool,
    num_threads: u32,
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("H", "", "set hashtable filename", "NAME");
    opts.optopt("t", "", "set taxonomy filename", "NAME");
    opts.optopt("o", "", "set options filename", "NAME");
    opts.optopt("O", "", "set output filename", "NAME");
    opts.optflag("m", "", "use MPA style output");
    opts.optflag("z", "", "report zeros");
    opts.optflag("s", "", "skip counts");
    opts.optopt("p", "", "set number of threads", "NUM");
    opts.optflag("h", "", "print this help menu");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f)
        }
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return;
    }
    let cmd_opts = CmdOptions {
        hashtable_filename: matches.opt_str("H").unwrap(),
        taxonomy_filename: matches.opt_str("t").unwrap(),
        options_filename: matches.opt_str("o").unwrap(),
        output_filename: matches.opt_str("O").unwrap_or("/dev/fd/1".to_string()),
        use_mpa_style: matches.opt_present("m"),
        report_zeros: matches.opt_present("z"),
        skip_counts: matches.opt_present("s"),
        num_threads: matches.opt_get_default("p", 1).unwrap(),
    };

    // Rest of your code here...
}
