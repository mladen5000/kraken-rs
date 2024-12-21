/*
* Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
*
* This file is part of the Kraken 2 taxonomic sequence classification system.
*/

#[cfg(not(feature = "openmp"))]
pub mod omp_hack {
    pub struct OmpLock;

    impl Default for OmpLock {
        fn default() -> Self {
            Self::new()
        }
    }

    impl OmpLock {
        pub fn new() -> Self {
            OmpLock
        }

        pub fn lock(&self) {}

        pub fn unlock(&self) {}

        pub fn try_lock(&self) -> bool {
            true
        }
    }

    pub fn get_thread_num() -> i32 {
        0
    }

    pub fn get_max_threads() -> i32 {
        1
    }

    pub fn set_num_threads(_num: i32) {}
}

#[cfg(feature = "openmp")]
pub mod omp_hack {
    use std::sync::Mutex;

    pub struct OmpLock(Mutex<()>);

    impl OmpLock {
        pub fn new() -> Self {
            OmpLock(Mutex::new(()))
        }

        pub fn lock(&self) {
            self.0.lock().unwrap();
        }

        pub fn unlock(&self) {
            self.0.unlock();
        }

        pub fn try_lock(&self) -> bool {
            self.0.try_lock().is_ok()
        }
    }

    pub fn get_thread_num() -> i32 {
        rayon::current_thread_index().try_into().unwrap()
    }

    pub fn get_max_threads() -> i32 {
        rayon::current_num_threads().try_into().unwrap()
    }

    pub fn set_num_threads(num: i32) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num as usize)
            .build_global()
            .expect("Failed to build global thread pool");
    }
}
