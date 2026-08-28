use engeom::Result;
use std::env;
use std::path::{Path, PathBuf};

/// The directory name searched for when locating the private test data, and the sibling directory
/// that test output is written into.
///
/// These carry the `engeom-` prefix deliberately. `find_private_test_data` walks up the parent
/// chain from the current working directory looking for the first of these, so a generic name like
/// `test-data` risks matching an unrelated directory that happens to sit in a shared parent.
const TEST_DATA_FOLDER: &str = "engeom-test-data";
const RESULT_DATA_FOLDER: &str = "engeom-result-data";

pub struct PathPair {
    root_path: PathBuf,
    suffix: Vec<String>,
}

impl PathPair {
    pub fn new_joined(&self, name: &str) -> Result<Self> {
        let mut items = self.suffix.clone();
        items.push(name.to_string());
        Self::try_from(&self.root_path, &items)
    }

    pub fn data_root(&self) -> PathBuf {
        self.root_path.join(TEST_DATA_FOLDER)
    }

    pub fn data(&self) -> PathBuf {
        let mut path = self.data_root();
        for s in self.suffix.iter() {
            path.push(s);
        }
        path
    }

    pub fn result_root(&self) -> PathBuf {
        self.root_path.join(RESULT_DATA_FOLDER)
    }

    pub fn result(&self) -> PathBuf {
        let mut path = self.result_root();
        for s in self.suffix.iter() {
            path.push(s);
        }
        path
    }

    pub fn try_from(root_path: &Path, items: &[String]) -> Result<PathPair> {
        let item = PathPair {
            root_path: root_path.to_path_buf(),
            suffix: items.to_vec(),
        };

        if !item.data().exists() {
            return Err(format!("Data path does not exist: {}", item.data().display()).into());
        }
        if !item.result().exists() {
            std::fs::create_dir_all(item.result())?;
        }

        Ok(item)
    }

    pub fn try_from_root(root_path: &Path) -> Result<PathPair> {
        Self::try_from(root_path, &[])
    }
}

pub fn find_private_test_data() -> Result<PathPair> {
    let mut current = env::current_dir()?;

    loop {
        let target = current.join(TEST_DATA_FOLDER);
        if target.exists() && target.is_dir() {
            return PathPair::try_from_root(&current);
        }

        current = current
            .parent()
            .ok_or("Could no longer recursively check parents")?
            .to_path_buf();
    }
}
