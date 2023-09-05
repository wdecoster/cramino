use rust_htslib::{bam, bam::Read};

pub fn get_genome_size(bamp: &String) -> u64 {
    let bam = bam::Reader::from_path(bamp).unwrap();
    let header = bam::Header::from_template(bam.header());
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
    genome_size
}

pub fn accuracy_to_phred(identity: f64) -> usize {
    // convert identity to phred scale
    // but return as usize (as that will be used for the histogram)
    // this is therefore not accurate for other applications
    (-10.0 * (1.0 - identity / 100.0).log10()) as usize
}
