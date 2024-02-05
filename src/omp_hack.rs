// Rust doesn't have a direct equivalent to OpenMP, but you can use the
// `rayon` crate for parallelism. If `rayon` isn't available, you can
// use the following fallback code:

pub struct OmpLock;

pub fn omp_get_thread_num() -> i32 {
    0
}

pub fn omp_get_max_threads() -> i32 {
    1
}

pub fn omp_set_num_threads(_num: i32) {}

pub fn omp_init_lock(_lock: &mut OmpLock) {}

pub fn omp_destroy_lock(_lock: &mut OmpLock) {}

pub fn omp_set_lock(_lock: &mut OmpLock) {}

pub fn omp_unset_lock(_lock: &mut OmpLock) {}

pub fn omp_test_lock(_lock: &mut OmpLock) -> i32 {
    1
}
