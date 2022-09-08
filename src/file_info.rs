use checksums::hash_file;
use chrono::offset::Local;
use chrono::DateTime;
use std::fmt;
use std::fs;
use std::path::Path;

pub(crate) struct BamFile {
    pub(crate) path: String,
}

impl BamFile {
    pub fn file_name(&self) -> String {
        Path::new(&self.path)
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string()
    }
    pub fn checksum(&self) -> String {
        hash_file(Path::new(&self.path), checksums::Algorithm::MD5)
    }

    pub fn file_time(&self) -> String {
        let metadata = fs::metadata(&self.path);
        if let Ok(time) = metadata.unwrap().created() {
            let datetime: DateTime<Local> = time.into();
            format!("{}", datetime.format("%d/%m/%Y %T"))
        } else {
            "NA".to_string()
        }
    }
}

impl fmt::Display for BamFile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.path)
    }
}
