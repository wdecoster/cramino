use std::path::PathBuf;

pub fn get_genome_size(
    header: &rust_htslib::bam::Header,
) -> Result<u64, rust_htslib::errors::Error> {
    let mut genome_size = 0;
    // print header records to the terminal, akin to samtool
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
pub fn calculate_data_yield(lengths: &[u128]) -> (u128, u128) {
    lengths.iter().fold((0u128, 0u128), |(total, long), &len| {
        let long_increment = if len > 25000 { len } else { 0 };
        (total + len, long + long_increment)
    })
}

pub fn is_file(pathname: &str) -> Result<(), String> {
    if pathname == "-"
        || pathname.starts_with("http")
        || pathname.starts_with("ftp")
        || pathname.starts_with("s3")
    {
        return Ok(());
    }
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}