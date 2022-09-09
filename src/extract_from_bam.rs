use log::error;
//use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::htslib;
use rust_htslib::{bam, bam::Read}; // for BAM_F*
use std::path::PathBuf;

pub fn extract(bamp: &String, threads: usize) -> (Vec<u32>, Vec<f32>) {
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(threads)
    //     .build()
    //     .unwrap();
    let mut bam = check_bam(bamp);
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    bam.fetch(".").expect("Failure fetching reads");
    bam.rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
        .map(|read| (read.seq_len() as u32, pid_from_cigar(read)))
        .unzip()
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

fn pid_from_cigar(r: std::rc::Rc<rust_htslib::bam::Record>) -> f32 {
    let mut matches = 0;
    let mut dels = 0;
    let mut ins = 0;
    for entry in r.cigar().iter() {
        match entry {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                matches += *len;
            }
            Cigar::Del(len) => {
                dels += *len;
            }
            Cigar::Ins(len) => {
                ins += *len;
            }
            _ => (),
        }
    }
    let alignment_length = matches + dels + ins;

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
