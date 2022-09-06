use checksums::hash_file;
use chrono::offset::Local;
use chrono::DateTime;
use std::fs;
use std::path::Path;

pub fn file_name(bamp: &String) -> String {
    let bampath = Path::new(&bamp);
    bampath.file_name().unwrap().to_str().unwrap().to_string()
}

pub fn checksum(bamp: &String) -> String {
    let bampath = Path::new(&bamp);
    hash_file(bampath, checksums::Algorithm::MD5)
}

pub fn file_time(bamp: &String) -> String {
    let metadata = fs::metadata(bamp);
    if let Ok(time) = metadata.unwrap().created() {
        let datetime: DateTime<Local> = time.into();
        format!("{}", datetime.format("%d/%m/%Y %T"))
    } else {
        "NA".to_string()
    }
}
