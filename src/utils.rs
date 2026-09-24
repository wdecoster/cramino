use std::path::PathBuf;

pub fn get_genome_size(
    header: &rust_htslib::bam::Header,
) -> Result<u64, rust_htslib::errors::Error> {
    let mut genome_size = 0;
    // sum of the lengths of the reference sequences (SQ lines) in the header
    for (key, records) in header.to_hashmap() {
        for record in records {
            if key == "SQ" {
                genome_size += record["LN"]
                    .parse::<u64>()
                    .expect("Failed parsing length of chromosomes");
            }
        }
    }
    Ok(genome_size)
}

pub fn accuracy_to_phred(identity: f64) -> usize {
    // convert identity to phred scale
    // but return as usize (as that will be used for the histogram)
    // this is therefore not accurate for other applications
    (-10.0 * (1.0 - identity / 100.0).log10()) as usize
}

// Helper function to calculate data yield
pub fn calculate_data_yield(lengths: &[u32]) -> (u128, u128) {
    lengths
        .iter()
        .map(|len| *len as u128)
        .fold((0u128, 0u128), |(total, long), len| {
            let long_increment = if len > 25000 { len } else { 0 };
            (total + len, long + long_increment)
        })
}

/// Remote input, which is opened as a URL
pub fn is_remote(input: &str) -> bool {
    ["s3://", "https://", "http://", "ftp://"]
        .iter()
        .any(|scheme| input.starts_with(scheme))
}

pub fn is_file(pathname: &str) -> Result<(), String> {
    if pathname == "-" || is_remote(pathname) {
        return Ok(());
    }
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!(
            "Input file {} does not exist or is not a file",
            path.display()
        ))
    }
}

/// Exits with a descriptive error message, for invalid input that cramino can't handle,
/// rather than panicking with a backtrace
pub fn exit_with_error(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1)
}
