use std::error::Error;

pub mod common;
pub mod func1;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn min_max(f0: f64, f1: f64) -> (f64, f64) {
    if f0 < f1 {
        (f0, f1)
    } else {
        (f1, f0)
    }
}

#[cfg(test)]
mod tests {}
