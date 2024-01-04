use std::error::Error;

pub mod common;
mod func1;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[cfg(test)]
mod tests {
    use super::*;
}
