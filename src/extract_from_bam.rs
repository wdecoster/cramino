use log::error;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::PathBuf;

pub fn extract(bamp: &String, every: usize) -> (Vec<u32>, Vec<u32>) {
    let mut bam = check_bam(bamp);

    bam.records()
        .map(|r| r.unwrap())
        .filter(check_read)
        .enumerate()
        .filter(|(index, _)| index % every == 0)
        .map(|(_, read)| extract_from_read(read))
        .unzip()
}

fn check_bam(bamp: &String) -> bam::IndexedReader {
    if !PathBuf::from(bamp).is_file() {
        error!("ERROR: path to bam file {} is not valid!\n\n", bamp);
        panic!();
    };
    match bam::IndexedReader::from_path(&bamp) {
        Ok(handle) => handle,
        Err(e) => {
            error!("Error opening BAM {}.\n{}", bamp, e);
            panic!();
        }
    }
}

fn check_read(read: &bam::Record) -> bool {
    !read.is_unmapped() && !read.is_secondary()
}

fn extract_from_read(read: bam::Record) -> (u32, u32) {
    (read.seq_len() as u32, pid_from_cigar(read))
}

fn pid_from_cigar(r: bam::Record) -> u32 {
    let mut lengths: HashMap<&str, u32> = HashMap::from([("match", 0), ("del", 0), ("ins", 0)]);
    for entry in r.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                *lengths.entry("match").or_insert(0) += *len;
            }
            Cigar::Del(len) => {
                *lengths.entry("del").or_insert(0) += *len;
            }
            Cigar::Ins(len) => {
                *lengths.entry("ins").or_insert(0) += *len;
            }
            _ => (),
        }
    }
    100 * (1 - get_nm_tag(&r)) / (lengths[&"match"] + lengths[&"del"] + lengths[&"ins"])
}

fn get_nm_tag(record: &bam::Record) -> u32 {
    match record.aux(b"NM") {
        Ok(value) => {
            if let Aux::U32(v) = value {
                v
            } else {
                panic!("Unexpected type of Aux {:?}", value)
            }
        }
        Err(_e) => panic!("Unexpected result while trying to access the NM tag"),
    }
}
