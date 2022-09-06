use log::error;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::PathBuf;

pub fn extract(bamp: &String, every: usize) -> (Vec<u32>, Vec<f32>) {
    let mut bam = check_bam(bamp);
    bam.fetch(".").unwrap();
    if every == 1 {
        bam.records()
            .map(|r| r.unwrap())
            .filter(check_read)
            .map(extract_from_read)
            .unzip()
    } else {
        bam.records()
            .enumerate()
            .filter(|(index, _)| index % every == 0)
            .map(|(_, r)| r.unwrap())
            .filter(check_read)
            .map(extract_from_read)
            .unzip()
    }
}

fn check_bam(bamp: &String) -> bam::IndexedReader {
    if !PathBuf::from(bamp).is_file() {
        error!("ERROR: path to bam file {bamp} is not valid!\n\n");
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

fn extract_from_read(read: bam::Record) -> (u32, f32) {
    (read.seq_len() as u32, pid_from_cigar(read))
}

fn pid_from_cigar(r: bam::Record) -> f32 {
    let mut operations: HashMap<&str, u32> = HashMap::from([("match", 0), ("del", 0), ("ins", 0)]);
    for entry in r.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                *operations.entry("match").or_insert(0) += *len;
            }
            Cigar::Del(len) => {
                *operations.entry("del").or_insert(0) += *len;
            }
            Cigar::Ins(len) => {
                *operations.entry("ins").or_insert(0) += *len;
            }
            _ => (),
        }
    }
    let alignment_length = operations[&"match"] + operations[&"del"] + operations[&"ins"];

    100.0 * (1.0 - get_nm_tag(&r) / alignment_length as f32)
}

fn get_nm_tag(record: &bam::Record) -> f32 {
    match record.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => v as f32,
            Aux::U16(v) => v as f32,
            Aux::U32(v) => v as f32,
            _ => panic!("Unexpected type of Aux {:?}", value),
        },
        Err(_e) => panic!("Unexpected result while trying to access the NM tag"),
    }
}
