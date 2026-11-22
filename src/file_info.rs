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
            .expect("Could not get file name")
            .to_str()
            .expect("Could not convert file name to string")
            .to_string()
    }

    pub fn file_time(&self) -> String {
        if self.path == "-" || crate::utils::is_remote(&self.path) {
            return "NA".to_string();
        }
        let metadata = fs::metadata(&self.path);
        if let Ok(time) = metadata
            .expect("Failed to extract metadata from file")
            .created()
        {
            match jiff::Zoned::try_from(time) {
                Ok(datetime) => datetime.strftime("%d/%m/%Y %T").to_string(),
                Err(_) => "NA".to_string(),
            }
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
