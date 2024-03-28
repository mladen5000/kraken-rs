use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};
use std::thread;

struct ThreadPool {
    workers: Vec<Worker>,
    sender: Sender<Message>,
}

struct Worker {
    id: usize,
    thread: Option<thread::JoinHandle<()>>,
}

enum Message {
    NewJob(Job),
    Terminate,
}

type Job = Box<dyn FnOnce() + Send + 'static>;

impl ThreadPool {
    fn new(size: usize) -> ThreadPool {
        assert!(size > 0);

        let (sender, receiver) = mpsc::channel();
        let receiver = Arc::new(Mutex::new(receiver));
        let mut workers = Vec::with_capacity(size);

        for id in 0..size {
            workers.push(Worker::new(id, Arc::clone(&receiver)));
        }

        ThreadPool { workers, sender }
    }

    fn execute<F>(&self, f: F)
    where
        F: FnOnce() + Send + 'static,
    {
        let job = Box::new(f);
        self.sender.send(Message::NewJob(job)).unwrap();
    }
}

impl Worker {
    fn new(id: usize, receiver: Arc<Mutex<Receiver<Message>>>) -> Worker {
        let thread = thread::spawn(move || loop {
            let message = receiver.lock().unwrap().recv().unwrap();
            match message {
                Message::NewJob(job) => {
                    job();
                }
                Message::Terminate => {
                    break;
                }
            }
        });

        Worker {
            id,
            thread: Some(thread),
        }
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        for _ in &self.workers {
            self.sender.send(Message::Terminate).unwrap();
        }

        for worker in &mut self.workers {
            if let Some(thread) = worker.thread.take() {
                thread.join().unwrap();
            }
        }
    }
}

fn main() {
    let pool = ThreadPool::new(4);
    pool.execute(|| {
        println!("Hello from the thread pool");
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};
    use std::time::Duration;

    #[test]
    fn test_thread_pool_execution() {
        let pool = ThreadPool::new(4);
        let data = Arc::new(Mutex::new(0));

        for _ in 0..8 {
            let data = Arc::clone(&data);
            pool.execute(move || {
                let mut num = data.lock().unwrap();
                *num += 1;
            });
        }

        // Allow some time for threads to execute jobs
        thread::sleep(Duration::from_secs(1));
        assert_eq!(*data.lock().unwrap(), 8);
    }

    #[test]
    fn test_thread_pool_multiple_executions() {
        let pool = ThreadPool::new(4);
        let data = Arc::new(Mutex::new(0));

        for _ in 0..4 {
            let data = Arc::clone(&data);
            pool.execute(move || {
                let mut num = data.lock().unwrap();
                *num += 2;
            });
        }

        // Allow some time for threads to execute jobs
        thread::sleep(Duration::from_secs(1));
        assert_eq!(*data.lock().unwrap(), 8);
    }

    #[test]
    #[should_panic(expected = "assertion failed")]
    fn test_thread_pool_zero_size() {
        ThreadPool::new(0); // Should panic due to zero size
    }
}
