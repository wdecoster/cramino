use bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};

pub struct Data {
    pub lengths: Option<Vec<u64>>,
    pub identities: Option<Vec<f64>>,
    pub tids: Option<Vec<i32>>,
    pub starts: Option<Vec<i64>>,
    pub ends: Option<Vec<i64>>,
    pub phasesets: Option<Vec<Option<u32>>>,
    pub exons: Option<Vec<usize>>,
}

pub fn extract(
    bam_path: &String,
    threads: usize,
    reference: Option<String>,
    min_read_len: usize,
    arrow: Option<String>,
    chroms: bool,
    phase: bool,
    spliced: bool,
) -> Data {
    let mut lengths = vec![];
    let mut identities = vec![];
    let mut tids = vec![];
    let mut starts = vec![];
    let mut ends = vec![];
    let mut phasesets = vec![];
    let mut exons = vec![];
    let mut bam = if bam_path == "-" {
        bam::Reader::from_stdin().expect("\n\nError reading alignments from stdin.\nDid you include the file header with -h?\n\n\n\n")
    } else {
        bam::Reader::from_path(bam_path)
            .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n")
    };
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");

    match reference {
        None => (),
        Some(s) => bam
            .set_reference(s)
            .expect("Failure setting bam/cram reference"),
    }
    for read in bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
        .filter(|read| read.seq_len() > min_read_len)
    {
        lengths.push(read.seq_len() as u64);
        if chroms || phase {
            tids.push(read.tid());
        }
        if phase {
            starts.push(read.pos());
            ends.push(read.reference_end());
            phasesets.push(get_phaseset(&read));
        }
        if spliced {
            exons.push(get_exon_number(&read));
        }
        identities.push(gap_compressed_identity(read));
    }
    match arrow {
        None => (),
        Some(s) => crate::feather::save_as_arrow(s, lengths.clone(), identities.clone()),
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    Data {
        lengths: Some(lengths),
        identities: Some(identities),
        tids: if chroms || phase { Some(tids) } else { None },
        starts: if phase { Some(starts) } else { None },
        ends: if phase { Some(ends) } else { None },
        phasesets: if phase { Some(phasesets) } else { None },
        exons: if spliced { Some(exons) } else { None },
    }
}

/// Calculates the gap-compressed identity
/// based on https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
/// recent minimap2 version have that as the de tag
/// if that is not present it is calculated from CIGAR and NM
fn gap_compressed_identity(record: std::rc::Rc<rust_htslib::bam::Record>) -> f64 {
    match get_de_tag(&record) {
        Some(v) => v as f64,
        None => {
            let mut matches = 0;
            let mut gap_size = 0;
            let mut gap_count = 0;
            for entry in record.cigar().iter() {
                match entry {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        matches += *len;
                    }
                    Cigar::Del(len) | Cigar::Ins(len) => {
                        gap_size += *len;
                        gap_count += 1;
                    }
                    _ => (),
                }
            }
            1.0 - ((get_nm_tag(&record) - gap_size + gap_count) as f64
                / (matches + gap_count) as f64)
        }
    }
}

fn get_nm_tag(record: &bam::Record) -> u32 {
    match record.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => u32::from(v),
            Aux::U16(v) => u32::from(v),
            Aux::U32(v) => v,
            Aux::I32(v) => u32::try_from(v).expect("Identified a negative NM tag"),
            _ => panic!("Unexpected type of Aux for NM tag: {:?}", value),
        },
        Err(_e) => panic!("Unexpected result while trying to access the NM tag"),
    }
}

/// Get the de:f tag from minimap2, which is the gap compressed sequence divergence
/// Which is converted into percent identity with 100 * (1 - de)
/// This tag can be absent if the aligner version is not quite recent
fn get_de_tag(record: &bam::Record) -> Option<f32> {
    match record.aux(b"de") {
        Ok(value) => match value {
            Aux::Float(v) => Some(100.0 * (1.0 - v)),
            _ => panic!("Unexpected type of Aux for de tag: {:?}", value),
        },
        Err(_e) => None,
    }
}

fn get_phaseset(record: &bam::Record) -> Option<u32> {
    match record.aux(b"PS") {
        Ok(value) => match value {
            Aux::U8(v) => Some(u32::from(v)),
            Aux::U16(v) => Some(u32::from(v)),
            Aux::U32(v) => Some(v),
            _ => panic!("Unexpected type of Aux for phaseset: {:?}", value),
        },
        Err(_e) => None,
    }
}

fn get_exon_number(record: &bam::Record) -> usize {
    let mut exon_count = 1;

    for op in record.cigar().iter() {
        if let Cigar::RefSkip(_len) = op {
            exon_count += 1;
        }
    }

    exon_count
}
