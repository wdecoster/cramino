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
}

// pub fn extract_unified(
//     bam_path: &String,
//     threads: usize,
//     min_read_len: usize,
//     arrow: Option<String>,
//     chroms: bool,
//     phase: bool,
// ) -> Data {
//     let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
//     bam.set_threads(threads)
//         .expect("Failure setting decompression threads");
//     use unzip_n::unzip_n;
//     unzip_n!(6);
//     let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
//     bam.set_threads(threads)
//         .expect("Failure setting decompression threads");
//     let (mut lengths, tids, starts, ends, phasesets, mut identities) = bam
//         .rc_records()
//         .map(|r| r.expect("Failure parsing Bam file"))
//         .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
//         .filter(|read| read.seq_len() > min_read_len)
//         .map(|read| {
//             (
//                 read.seq_len() as u64,
//                 if chroms { Some(read.tid()) } else { None },
//                 if phase { Some(read.pos()) } else { None },
//                 read.reference_end(),
//                 get_phaseset(&read),
//                 gap_compressed_identity(read),
//             )
//         })
//         .unzip_n_vec();
//     match arrow {
//         None => (),
//         Some(s) => crate::feather::save_as_arrow(s, lengths.clone(), identities.clone()),
//     }
//     lengths.sort_unstable();
//     identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
//     Data {
//         lengths: Some(lengths),
//         identities: Some(identities),
//         tids: Some(tids),
//         starts: Some(starts),
//         ends: Some(ends),
//         phasesets: Some(phasesets),
//     }
// }

pub fn extract(
    bam_path: &String,
    threads: usize,
    min_read_len: usize,
    arrow: Option<String>,
) -> Data {
    let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    let (mut lengths, mut identities): (Vec<u64>, Vec<f64>) = bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
        .filter(|read| read.seq_len() > min_read_len)
        .map(|read| (read.seq_len() as u64, gap_compressed_identity(read)))
        .unzip();
    match arrow {
        None => (),
        Some(s) => crate::feather::save_as_arrow(s, lengths.clone(), identities.clone()),
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    Data {
        lengths: Some(lengths),
        identities: Some(identities),
        tids: None,
        starts: None,
        ends: None,
        phasesets: None,
    }
}

pub fn extract_with_chroms(
    bam_path: &String,
    threads: usize,
    min_read_len: usize,
    arrow: Option<String>,
) -> Data {
    use unzip_n::unzip_n;
    unzip_n!(3);
    let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    let (mut lengths, tids, mut identities) = bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
        .filter(|read| read.seq_len() > min_read_len)
        .map(|read| {
            (
                read.seq_len() as u64,
                read.tid(),
                gap_compressed_identity(read),
            )
        })
        .unzip_n_vec();
    match arrow {
        None => (),
        Some(s) => crate::feather::save_as_arrow(s, lengths.clone(), identities.clone()),
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    Data {
        lengths: Some(lengths),
        identities: Some(identities),
        tids: Some(tids),
        starts: None,
        ends: None,
        phasesets: None,
    }
}

pub fn extract_with_phase(
    bam_path: &String,
    threads: usize,
    min_read_len: usize,
    arrow: Option<String>,
) -> Data {
    use unzip_n::unzip_n;
    unzip_n!(6);
    let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    let (mut lengths, tids, starts, ends, phasesets, mut identities) = bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
        .filter(|read| read.seq_len() > min_read_len)
        .map(|read| {
            (
                read.seq_len() as u64,
                read.tid(),
                read.pos(),
                read.reference_end(),
                get_phaseset(&read),
                gap_compressed_identity(read),
            )
        })
        .unzip_n_vec();
    match arrow {
        None => (),
        Some(s) => crate::feather::save_as_arrow(s, lengths.clone(), identities.clone()),
    }
    lengths.sort_unstable();
    identities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    Data {
        lengths: Some(lengths),
        identities: Some(identities),
        tids: Some(tids),
        starts: Some(starts),
        ends: Some(ends),
        phasesets: Some(phasesets),
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
            _ => panic!("Unexpected type of Aux {:?}", value),
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
            _ => panic!("Unexpected type of Aux {:?}", value),
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
            _ => panic!("Unexpected type of Aux {:?}", value),
        },
        Err(_e) => None,
    }
}
