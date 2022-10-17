use rust_htslib::{bam, bam::Read};

pub fn get_genome_size(bamp: &String) -> u64 {
    let bam = bam::Reader::from_path(bamp).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut genome_size = 0;
    // print header records to the terminal, akin to samtool
    for (key, records) in header.to_hashmap() {
        for record in records {
            if key == "SQ" {
            genome_size += record["LN"].parse::<u64>().expect("Failed parsing length of chromosomes");
        }
    }
    }
    genome_size
}
