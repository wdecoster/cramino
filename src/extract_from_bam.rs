use bam::ext::BamRecordExtensions;
use log::warn;
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};
use std::env;
use url::Url;

pub struct Data {
    /// Aligned length (without clipped bases) of all primary and supplementary alignments,
    /// and the full length of unmapped reads with --ubam
    pub lengths: Option<Vec<u32>>,
    /// Full length of each read, i.e. only of the primary alignments (and unmapped reads with --ubam)
    pub read_lengths: Vec<u32>,
    pub num_reads: usize,
    /// Number of reads in the file, before filtering: primary alignments and unmapped reads
    pub all_reads: usize,
    pub identities: Option<Vec<f64>>,
    pub q_score_hist: Option<QScoreHistogramData>,
    /// Chromosome of each alignment, including supplementary alignments, for the phaseblocks
    pub tids: Option<Vec<i32>>,
    /// Chromosome of each read (primary alignment), for the karyotype
    pub karyotype_tids: Option<Vec<i32>>,
    pub starts: Option<Vec<i64>>,
    pub ends: Option<Vec<i64>>,
    /// Phaseset of each alignment, including supplementary alignments, to determine the phaseblocks
    pub phasesets: Option<Vec<Option<u32>>>,
    /// Number of reads (primary alignments) with a phaseset
    pub num_phased_reads: usize,
    pub exons: Option<Vec<usize>>,
    pub is_ubam: bool,
}

pub struct QScoreHistogramData {
    pub counts: Vec<u64>,
    pub bases: Vec<u128>,
}

/// Sets up the CURL_CA_BUNDLE environment variable for HTTPS/S3 access
/// Tries to use a CA bundle from standard locations, with appropriate fallbacks
///
/// # Safety
/// Modifying the environment is only sound while no other threads can read it,
/// so this has to be called at the start of main, before any threads are spawned
/// (by htslib, rayon or anything else).
pub unsafe fn setup_ssl_certificates() {
    // Only configure if not already set by the user
    if env::var("CURL_CA_BUNDLE").is_ok() {
        return;
    }

    // Common CA bundle locations across different systems
    let possible_paths = vec![
        "/etc/ssl/certs/ca-certificates.crt",     // Debian/Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",       // RHEL/CentOS/Amazon Linux
        "/etc/ssl/ca-bundle.pem",                 // SUSE
        "/usr/local/share/certs/ca-root-nss.crt", // FreeBSD
        "/usr/local/etc/openssl/cert.pem",        // macOS Homebrew
        "/etc/ssl/cert.pem",                      // macOS/OpenBSD
    ];

    // Try each path in order
    for path in possible_paths {
        if std::path::Path::new(path).exists() {
            // SAFETY: the caller guarantees that the program is still single-threaded
            unsafe { env::set_var("CURL_CA_BUNDLE", path) };
            return;
        }
    }

    // None of the paths exist, warn the user
    warn!(
        "Could not find a valid CA certificate bundle for HTTPS/S3 access. \
        HTTPS/S3 connections may fail. Set the CURL_CA_BUNDLE environment \
        variable to the path of your system's CA certificate bundle."
    );
}

pub fn extract(args: &crate::Cli) -> (Data, rust_htslib::bam::Header) {
    let mut lengths = vec![];
    let mut read_lengths = vec![];
    let mut num_reads = 0;
    let mut identities = vec![];
    let hist_requested = args.hist.is_some() || args.hist_count.is_some();
    let mut q_score_counts = Vec::new();
    let mut q_score_bases = Vec::new();
    if hist_requested {
        q_score_counts = vec![0u64; 41];
        q_score_bases = vec![0u128; 41];
    }
    let mut tids = vec![];
    let mut karyotype_tids = vec![];
    let mut starts = vec![];
    let mut ends = vec![];
    let mut phasesets = vec![];
    let mut num_phased_reads = 0;
    let mut exons = vec![];
    use crate::utils::exit_with_error;
    let mut bam = if args.input == "-" {
        bam::Reader::from_stdin().unwrap_or_else(|err| {
            exit_with_error(&format!(
                "Could not read alignments from stdin ({err}).\nDid you include the header, e.g. with samtools view -h?"
            ))
        })
    } else if crate::utils::is_remote(&args.input) {
        // the certificates for remote access are set up at the start of main
        let url = Url::parse(&args.input)
            .unwrap_or_else(|err| exit_with_error(&format!("Invalid URL {} ({err})", args.input)));
        bam::Reader::from_url(&url).unwrap_or_else(|err| exit_with_error(&err.to_string()))
    } else {
        bam::Reader::from_path(&args.input)
            .unwrap_or_else(|err| exit_with_error(&format!("{err}, is it a valid BAM/CRAM file?")))
    };
    // the reference is also set for other input than a .cram file, e.g. a CRAM on stdin.
    // It is not used for BAM input
    if let Some(reference) = &args.reference {
        bam.set_reference(reference).unwrap_or_else(|err| {
            exit_with_error(&format!("Could not set reference {reference} ({err})"))
        });
    }
    if args.input.ends_with(".cram") {
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_CIGAR
                | hts_sys::sam_fields_SAM_SEQ,
        )
        .expect("Failed setting cram options");
    }
    let header = bam.header().clone();
    let header = rust_htslib::bam::Header::from_template(&header);
    bam.set_threads(args.threads)
        .expect("Failure setting decompression threads");

    let min_read_len = u32::try_from(args.min_read_len).unwrap_or(u32::MAX);
    let mut all_reads = 0;
    for read in bam
        .rc_records()
        .map(|r| {
            r.unwrap_or_else(|err| {
                exit_with_error(&format!(
                    "Could not read a record from {} ({err}). For CRAM, this can be caused by a missing reference, which can be provided with --reference.",
                    args.input
                ))
            })
        })
        .filter(|record| record.flags() & (htslib::BAM_FSECONDARY) as u16 == 0)
        // every read has exactly one record that is neither secondary nor supplementary
        .inspect(|record| {
            if !record.is_supplementary() {
                all_reads += 1
            }
        })
        // unmapped reads are only kept with --ubam
        .filter(|record| args.ubam || record.flags() & htslib::BAM_FUNMAP as u16 == 0)
    {
        let read_length = cigar_query_length(&read, ALIGNED_OPS);
        if read_length < min_read_len {
            continue;
        }
        lengths.push(read_length);
        if !read.is_supplementary() {
            num_reads += 1;
            read_lengths.push(cigar_query_length(&read, READ_OPS));
            if args.karyotype {
                karyotype_tids.push(read.tid());
            }
        }
        if args.phased {
            tids.push(read.tid());
        }
        if args.phased {
            starts.push(read.pos());
            ends.push(read.reference_end());
            let phaseset = get_phaseset(&read);
            if phaseset.is_some() && !read.is_supplementary() {
                num_phased_reads += 1;
            }
            phasesets.push(phaseset);
        }
        if args.spliced {
            exons.push(get_exon_number(&read));
        }
        // with --ubam, the identity is estimated from the base qualities
        let identity = if args.ubam {
            qscore_to_accuracy(&read)
        } else {
            gap_compressed_identity(&read)
        };
        if let Some(identity) = identity
            && hist_requested
        {
            let phred = crate::utils::accuracy_to_phred(identity);
            let index = if phred < 40 { phred } else { 40 };
            q_score_counts[index] += 1;
            q_score_bases[index] += read_length as u128;
        }
        // reads without an identity are kept as NaN, so that the identities stay in the
        // same order as the lengths for the arrow output, and are left out afterwards
        identities.push(identity.unwrap_or(f64::NAN));
    }
    if let Some(s) = &args.arrow {
        crate::feather::save_as_arrow(
            s.to_string(),
            lengths.iter().map(|x| *x as u64).collect(),
            &identities,
        );
    }

    // sort vectors in descending order (required for N50/N75)
    lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    read_lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    identities.retain(|identity| !identity.is_nan());
    identities.par_sort_unstable_by(|a, b| b.total_cmp(a));
    (
        Data {
            lengths: Some(lengths),
            read_lengths,
            num_reads,
            all_reads,
            identities: Some(identities),
            q_score_hist: if hist_requested {
                Some(QScoreHistogramData {
                    counts: q_score_counts,
                    bases: q_score_bases,
                })
            } else {
                None
            },
            tids: if args.phased { Some(tids) } else { None },
            karyotype_tids: if args.karyotype {
                Some(karyotype_tids)
            } else {
                None
            },
            starts: if args.phased { Some(starts) } else { None },
            ends: if args.phased { Some(ends) } else { None },
            phasesets: if args.phased { Some(phasesets) } else { None },
            num_phased_reads,
            exons: if args.spliced { Some(exons) } else { None },
            is_ubam: args.ubam,
        },
        header,
    )
}

/// Calculates the gap-compressed identity
/// based on https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
/// recent minimap2 version have that as the de tag
/// if that is not present it is calculated from CIGAR and NM
/// None if the identity can't be determined, i.e. without aligned bases
fn gap_compressed_identity(record: &bam::Record) -> Option<f64> {
    let identity = match get_de_tag(record) {
        Some(v) => v,
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
            // NM includes the inserted and deleted bases, so should be at least gap_size.
            // An inconsistent, lower NM is taken as no mismatches rather than underflowing
            let mismatches = get_nm_tag(record).saturating_sub(gap_size);
            if matches + gap_count == 0 {
                return None;
            }
            100.0 * (1.0 - ((mismatches + gap_count) as f64 / (matches + gap_count) as f64))
        }
    };
    identity.is_finite().then_some(identity)
}

/// Computes the mean estimated accuracy (%) from per-base Q-scores
/// Q-score is Phred-scaled: Q = -10 * log10(P_error)
/// This function converts each Q-score to probability of correctness
/// and returns the average as a percentage.
/// None for reads without base qualities, which are therefore left out
fn qscore_to_accuracy(record: &bam::Record) -> Option<f64> {
    let quals = record.qual();
    if quals.is_empty() || quals.iter().all(|&q| q == 255) {
        // 255 indicates missing quality
        return None;
    }

    let sum_accuracy: f64 = quals
        .iter()
        .map(|&q| {
            // P_error = 10^(-Q/10), P_correct = 1 - P_error
            1.0 - 10_f64.powf(-(q as f64) / 10.0)
        })
        .sum();

    Some(100.0 * sum_accuracy / quals.len() as f64)
}

/// The NM tag is optional according to the SAM specification, but required here to calculate
/// the identity if the de tag is absent. Without it, cramino exits with an error.
fn get_nm_tag(record: &bam::Record) -> u32 {
    let nm = match record.aux(b"NM") {
        Ok(Aux::U8(v)) => Some(u32::from(v)),
        Ok(Aux::U16(v)) => Some(u32::from(v)),
        Ok(Aux::U32(v)) => Some(v),
        Ok(Aux::I8(v)) => u32::try_from(v).ok(),
        Ok(Aux::I16(v)) => u32::try_from(v).ok(),
        Ok(Aux::I32(v)) => u32::try_from(v).ok(),
        Ok(value) => crate::utils::exit_with_error(&format!(
            "Read {} has an NM tag of an unexpected type ({value:?}), while an integer is required.",
            String::from_utf8_lossy(record.qname())
        )),
        Err(_) => crate::utils::exit_with_error(&format!(
            "Read {} has neither an NM nor a de tag, one of which is required to calculate the identity of aligned reads.\n\
            NM tags can be added with `samtools calmd`.",
            String::from_utf8_lossy(record.qname())
        )),
    };
    nm.unwrap_or_else(|| {
        crate::utils::exit_with_error(&format!(
            "Read {} has a negative NM tag.",
            String::from_utf8_lossy(record.qname())
        ))
    })
}

/// Get the de:f tag from minimap2, which is the gap compressed sequence divergence
/// Which is converted into percent identity with 100 * (1 - de)
/// This tag can be absent if the aligner version is not quite recent.
/// A de tag of another type than float or double is ignored, falling back to the NM tag
fn get_de_tag(record: &bam::Record) -> Option<f64> {
    match record.aux(b"de") {
        Ok(Aux::Float(v)) => Some(100.0 * (1.0 - v as f64)),
        Ok(Aux::Double(v)) => Some(100.0 * (1.0 - v)),
        _ => None,
    }
}

/// The phaseset from the PS tag, which has to be a non-negative integer.
/// Otherwise cramino exits with an error
fn get_phaseset(record: &bam::Record) -> Option<u32> {
    let phaseset = match record.aux(b"PS") {
        Err(_) => return None,
        Ok(Aux::U8(v)) => Some(u32::from(v)),
        Ok(Aux::U16(v)) => Some(u32::from(v)),
        Ok(Aux::U32(v)) => Some(v),
        Ok(Aux::I8(v)) => u32::try_from(v).ok(),
        Ok(Aux::I16(v)) => u32::try_from(v).ok(),
        Ok(Aux::I32(v)) => u32::try_from(v).ok(),
        Ok(value) => crate::utils::exit_with_error(&format!(
            "Read {} has a PS tag of an unexpected type ({value:?}), while an integer is required.",
            String::from_utf8_lossy(record.qname())
        )),
    };
    Some(phaseset.unwrap_or_else(|| {
        crate::utils::exit_with_error(&format!(
            "Read {} has a negative PS tag.",
            String::from_utf8_lossy(record.qname())
        ))
    }))
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

/// CIGAR operations for the query bases in the alignment
const ALIGNED_OPS: u32 = 1 << htslib::BAM_CMATCH
    | 1 << htslib::BAM_CINS
    | 1 << htslib::BAM_CEQUAL
    | 1 << htslib::BAM_CDIFF;
/// CIGAR operations for all bases of the read, including those clipped from the alignment
const READ_OPS: u32 = ALIGNED_OPS | 1 << htslib::BAM_CSOFT_CLIP | 1 << htslib::BAM_CHARD_CLIP;

/// Number of query bases in the CIGAR operations selected by the `ops` bitmask.
/// This is taken from the CIGAR rather than from the sequence, as the latter can be absent (SEQ '*'),
/// and hard clipped bases are not in the sequence either.
/// The packed CIGAR is read directly, as unpacking it would allocate for every record.
/// Records without a CIGAR (unmapped reads, or rarely mapped reads with CIGAR '*') have their full sequence length.
/// Lengths are stored as u32 to save memory, which fits reads up to 4.29 Gbp.
fn cigar_query_length(read: &bam::Record, ops: u32) -> u32 {
    let cigar = read.raw_cigar();
    let length = if cigar.is_empty() {
        read.seq_len() as u64
    } else {
        cigar
            .iter()
            .filter(|op| ops & (1 << (*op & htslib::BAM_CIGAR_MASK)) != 0)
            .map(|op| (op >> htslib::BAM_CIGAR_SHIFT) as u64)
            .sum()
    };
    u32::try_from(length).expect("Read length does not fit in 32 bits")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test the Q-score to accuracy conversion formula
    /// Q = 10 means P_error = 0.1, P_correct = 0.9, accuracy = 90%
    /// Q = 20 means P_error = 0.01, P_correct = 0.99, accuracy = 99%
    /// Q = 30 means P_error = 0.001, P_correct = 0.999, accuracy = 99.9%
    #[test]
    fn test_qscore_to_probability_formula() {
        // Q = 10: error prob = 10^(-10/10) = 0.1, accuracy = 0.9
        let q10_accuracy = 1.0 - 10_f64.powf(-10.0 / 10.0);
        assert!((q10_accuracy - 0.9).abs() < 1e-10);

        // Q = 20: error prob = 10^(-20/10) = 0.01, accuracy = 0.99
        let q20_accuracy = 1.0 - 10_f64.powf(-20.0 / 10.0);
        assert!((q20_accuracy - 0.99).abs() < 1e-10);

        // Q = 30: error prob = 10^(-30/10) = 0.001, accuracy = 0.999
        let q30_accuracy = 1.0 - 10_f64.powf(-30.0 / 10.0);
        assert!((q30_accuracy - 0.999).abs() < 1e-10);
    }

    #[test]
    fn test_qscore_to_accuracy_with_record() {
        // Create a test record with known quality scores
        let mut record = bam::Record::new();
        // Set a simple read with quality scores of Q20
        // For a read with all Q20 bases, expected accuracy = 99%
        let qname = b"test_read";
        let seq = b"ACGT";
        let qual = vec![20u8; 4]; // All Q20
        record.set(qname, None, seq, &qual);

        let accuracy = qscore_to_accuracy(&record).unwrap();
        // Expected: 100 * (1 - 0.01) = 99.0
        assert!((accuracy - 99.0).abs() < 0.01);
    }

    #[test]
    fn test_qscore_to_accuracy_mixed_qualities() {
        // Create a test record with mixed quality scores
        let mut record = bam::Record::new();
        let qname = b"test_read";
        let seq = b"ACGT";
        // Mix of Q10, Q20, Q30, Q40
        let qual = vec![10u8, 20, 30, 40];
        record.set(qname, None, seq, &qual);

        let accuracy = qscore_to_accuracy(&record).unwrap();
        // Q10: 0.9, Q20: 0.99, Q30: 0.999, Q40: 0.9999
        // Average: (0.9 + 0.99 + 0.999 + 0.9999) / 4 = 0.972225
        // As percentage: 97.2225
        let expected = 100.0 * (0.9 + 0.99 + 0.999 + 0.9999) / 4.0;
        assert!((accuracy - expected).abs() < 0.01);
    }

    #[test]
    fn test_qscore_to_accuracy_missing_quality() {
        // Create a test record with missing quality (all 255)
        let mut record = bam::Record::new();
        let qname = b"test_read";
        let seq = b"ACGT";
        let qual = vec![255u8; 4]; // Missing quality indicator
        record.set(qname, None, seq, &qual);

        // reads without base qualities have no identity, and are left out
        assert_eq!(qscore_to_accuracy(&record), None);
    }

    fn record_with_cigar(seq: &[u8], cigar: &str) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = bam::record::CigarString::try_from(cigar).unwrap();
        let qual = vec![20u8; seq.len()];
        record.set(b"test_read", Some(&cigar), seq, &qual);
        record
    }

    #[test]
    fn test_query_length_with_and_without_clips() {
        // soft clips (primary) and hard clips (supplementary) are only part of the read length,
        // deletions don't consume query bases
        let record = record_with_cigar(&[b'A'; 80], "10S50M5I5D15S");
        assert_eq!(cigar_query_length(&record, ALIGNED_OPS), 55);
        assert_eq!(cigar_query_length(&record, READ_OPS), 80);
        let record = record_with_cigar(&[b'A'; 55], "100H50M5I5D200H");
        assert_eq!(cigar_query_length(&record, ALIGNED_OPS), 55);
        assert_eq!(cigar_query_length(&record, READ_OPS), 355);
    }

    #[test]
    fn test_query_length_without_sequence() {
        // SEQ '*' with soft clips used to underflow, as the clips were subtracted from a zero sequence length
        let record = record_with_cigar(b"", "10S50M10S");
        assert_eq!(cigar_query_length(&record, ALIGNED_OPS), 50);
        assert_eq!(cigar_query_length(&record, READ_OPS), 70);
    }

    #[test]
    fn test_query_length_unmapped() {
        let mut record = bam::Record::new();
        record.set(b"test_read", None, &[b'A'; 30], &[20u8; 30]);
        assert_eq!(cigar_query_length(&record, ALIGNED_OPS), 30);
        assert_eq!(cigar_query_length(&record, READ_OPS), 30);
    }

    #[test]
    fn test_identity_with_nm_lower_than_indel_bases() {
        // NM:i:0 is inconsistent with a 100 bp deletion, and used to underflow
        let mut record = record_with_cigar(&[b'A'; 1000], "500M100D500M");
        record.push_aux(b"NM", Aux::U8(0)).unwrap();
        let identity = gap_compressed_identity(&record).unwrap();
        // the deletion counts as one gap: 1 - 1 / (1000 + 1)
        assert!((identity - 100.0 * (1.0 - 1.0 / 1001.0)).abs() < 1e-9);
    }

    #[test]
    fn test_identity_without_aligned_bases() {
        // a mapped record without CIGAR has no aligned bases, and therefore no identity
        let mut record = record_with_cigar(&[b'A'; 10], "");
        record.push_aux(b"NM", Aux::U8(0)).unwrap();
        assert_eq!(gap_compressed_identity(&record), None);
    }

    #[test]
    fn test_identity_from_de_tag_as_float_or_double() {
        let mut record = record_with_cigar(&[b'A'; 10], "10M");
        record.push_aux(b"de", Aux::Double(0.02)).unwrap();
        assert!((gap_compressed_identity(&record).unwrap() - 98.0).abs() < 1e-9);
        let mut record = record_with_cigar(&[b'A'; 10], "10M");
        record.push_aux(b"de", Aux::Float(0.02)).unwrap();
        assert!((gap_compressed_identity(&record).unwrap() - 98.0).abs() < 1e-6);
    }
}
