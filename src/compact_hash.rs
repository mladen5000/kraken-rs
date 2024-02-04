pub use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

pub struct CompactHashCell {
    data: u32,
}

impl CompactHashCell {
    fn hashed_key(&self, value_bits: usize) -> u64 {
        self.data as u64 >> value_bits
    }

    fn value(&self, value_bits: usize) -> u64 {
        self.data as u64 & ((1 << value_bits) - 1)
    }

    fn populate(&mut self, compacted_key: u64, val: u64, key_bits: usize, value_bits: usize) {
        assert_eq!(
            key_bits + value_bits,
            32,
            "key len of {} and value len of {} don't sum to 32",
            key_bits,
            value_bits
        );
        assert!(
            key_bits > 0 && value_bits > 0,
            "key len and value len must be nonzero"
        );
        let max_value = (1u64 << value_bits) - 1;
        assert!(
            max_value >= val,
            "value len of {} too small for value of {}",
            value_bits,
            val
        );
        self.data = ((compacted_key << value_bits) | val) as u32;
    }
}

pub struct CompactHashTable {
    capacity: AtomicUsize,
    size: AtomicUsize,
    key_bits: usize,
    value_bits: usize,
    table: Mutex<Vec<CompactHashCell>>,
    file_backed: bool,
    locks_initialized: bool,
    backing_file: String,
    zone_locks: [Mutex<u8>; 256],
}

pub trait CompactHashTableTrait {
    fn new(capacity: usize, key_bits: usize, value_bits: usize) -> Self;
    fn from_file(filename: &str, memory_mapping: bool) -> Self;
    fn get(&self, key: u64) -> Option<u64>;
    fn find_index(&self, key: u64) -> Option<usize>;
    fn compare_and_set(&self, key: u64, new_value: u64, old_value: &mut u64) -> bool;
    fn direct_compare_and_set(
        &self,
        idx: usize,
        key: u64,
        new_value: u64,
        old_value: &mut u64,
    ) -> bool;
    fn write_table(&self, filename: &str);
    fn get_value_counts(&self) -> HashMap<u64, u64>;
    fn capacity(&self) -> usize;
    fn size(&self) -> usize;
    fn key_bits(&self) -> usize;
    fn value_bits(&self) -> usize;
    fn occupancy(&self) -> f64;
    fn second_hash(&self, first_hash: u64) -> u64;
}

impl CompactHashTableTrait for CompactHashTable {
    // Methods go here
    fn new(capacity: usize, key_bits: usize, value_bits: usize) -> Self {
        assert_eq!(
            key_bits + value_bits,
            32,
            "sum of key bits and value bits must equal {}",
            32
        );
        assert!(key_bits > 0, "key bits cannot be zero");
        assert!(value_bits > 0, "value bits cannot be zero");

        let mut table = vec![CompactHashCell { data: 0 }; capacity];

        CompactHashTable {
            capacity: AtomicUsize::new(capacity),
            size: AtomicUsize::new(0),
            key_bits,
            value_bits,
            table: Mutex::new(table),
            file_backed: false,
            locks_initialized: true,
            backing_file: String::new(),
            zone_locks: Default::default(),
        }
    }
    fn from_file(filename: &str, memory_mapping: bool) -> Self {
        let mut table = CompactHashTable::new(0, 0, 0);
        table.load_table(filename, memory_mapping);
        table
    }
    fn get(&self, key: u64) -> Option<u64> {
        use crate::kv_store::murmurhash3;
        let hc = murmurhash3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity.load(Ordering::Relaxed) as u64;
        let first_idx = idx;
        let mut step = 0;
        let table = self.table.lock().unwrap();

        loop {
            if table[idx as usize].value(self.value_bits) == 0 {
                break;
            }
            if table[idx as usize].hashed_key(self.value_bits) == compacted_key {
                return Some(table[idx as usize].value(self.value_bits));
            }
            if step == 0 {
                step = second_hash(hc);
            }
            idx += step;
            idx %= self.capacity.load(Ordering::Relaxed) as u64;
            if idx == first_idx {
                break;
            }
        }
        None
    }
    fn find_index(&self, key: u64) -> Option<usize> {
        let hc = murmur_hash_3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity.load(Ordering::Relaxed) as u64;
        let first_idx = idx;
        let mut step = 0;
        let table = self.table.lock().unwrap();

        loop {
            if table[idx as usize].value(self.value_bits) == 0 {
                return None;
            }
            if table[idx as usize].hashed_key(self.value_bits) == compacted_key {
                return Some(idx as usize);
            }
            if step == 0 {
                step = second_hash(hc);
            }
            idx += step;
            idx %= self.capacity.load(Ordering::Relaxed) as u64;
            if idx == first_idx {
                break;
            }
        }
        None
    }
    fn compare_and_set(&self, key: u64, new_value: u64, old_value: &mut u64) -> bool {
        if self.file_backed || new_value == 0 {
            return false;
        }
        let hc = murmur_hash_3(key);
        let compacted_key = hc >> (32 + self.value_bits);
        let mut idx = hc % self.capacity.load(Ordering::Relaxed) as u64;
        let first_idx = idx;
        let mut step = 0;
        let mut set_successful = false;
        let mut search_successful = false;
        let mut table = self.table.lock().unwrap();

        while !search_successful {
            if table[idx as usize].value(self.value_bits) == 0
                || table[idx as usize].hashed_key(self.value_bits) == compacted_key
            {
                search_successful = true;
                if *old_value == table[idx as usize].value(self.value_bits) {
                    table[idx as usize].populate(
                        compacted_key,
                        new_value,
                        self.key_bits,
                        self.value_bits,
                    );
                    if *old_value == 0 {
                        self.size.fetch_add(1, Ordering::SeqCst);
                    }
                    set_successful = true;
                } else {
                    *old_value = table[idx as usize].value(self.value_bits);
                }
            }
            if step == 0 {
                step = second_hash(hc);
            }
            idx += step;
            idx %= self.capacity.load(Ordering::Relaxed) as u64;
            if idx == first_idx {
                panic!("compact hash table capacity exceeded");
            }
        }
        set_successful
    }
    fn direct_compare_and_set(
        &self,
        idx: usize,
        key: u64,
        new_value: u64,
        old_value: &mut u64,
    ) -> bool {
    }
    fn write_table(&self, filename: &str) {
        let mut ofs = std::fs::File::create(filename).expect("Unable to create file");
        let capacity = self.capacity.load(Ordering::Relaxed);
        let size = self.size.load(Ordering::Relaxed);
        let key_bits = self.key_bits;
        let value_bits = self.value_bits;
        let table = self.table.lock().unwrap();

        ofs.write_all(&capacity.to_ne_bytes())
            .expect("Unable to write data");
        ofs.write_all(&size.to_ne_bytes())
            .expect("Unable to write data");
        ofs.write_all(&key_bits.to_ne_bytes())
            .expect("Unable to write data");
        ofs.write_all(&value_bits.to_ne_bytes())
            .expect("Unable to write data");

        for cell in table.iter() {
            ofs.write_all(&cell.data.to_ne_bytes())
                .expect("Unable to write data");
        }
    }
    fn get_value_counts(&self) -> HashMap<u64, u64> {
        let table = self.table.lock().unwrap();
        let mut value_counts: HashMap<u64, u64> = HashMap::new();

        table.par_iter().for_each(|cell| {
            let val = cell.value(self.value_bits);
            if val != 0 {
                *value_counts.entry(val).or_insert(0) += 1;
            }
        });

        value_counts
    }
    fn capacity(&self) -> usize {
        self.capacity
    }
    fn size(&self) -> usize {
        self.size
    }
    fn key_bits(&self) -> usize {
        self.key_bits
    }
    fn value_bits(&self) -> usize {
        self.value_bits
    }
    fn occupancy(&self) -> f64 {
        self.size * 1.0 / self.capacity
    }
    fn second_hash(&self, first_hash: u64) -> u64 {
        #[cfg(feature = "linear_probing")]
        {
            1
        }
        #[cfg(not(feature = "linear_probing"))]
        {
            (first_hash >> 8) | 1
        }
    }
}
