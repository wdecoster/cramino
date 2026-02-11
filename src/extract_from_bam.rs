use bam::ext::BamRecordExtensions;
use log::warn;
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};
use std::env;
use url::Url;

pub struct Data {
    pub lengths: Option<Vec<u128>>,
    pub num_reads: usize,
    pub all_counts: usize,
    pub identities: Option<Vec<f64>>,
    pub q_score_hist: Option<QScoreHistogramData>,
    pub tids: Option<Vec<i32>>,
    pub starts: Option<Vec<i64>>,
    pub ends: Option<Vec<i64>>,
    pub phasesets: Option<Vec<Option<u32>>>,
    pub exons: Option<Vec<usize>>,
    pub is_ubam: bool,
}

pub struct QScoreHistogramData {
    pub counts: Vec<u64>,
    pub bases: Vec<u128>,
}

/// Sets up the CURL_CA_BUNDLE environment variable for HTTPS/S3 access
/// Tries to use a CA bundle from standard locations, with appropriate fallbacks
fn setup_ssl_certificates() {
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
            // TODO: Audit that the environment access only happens in single-threaded code.
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
    let mut starts = vec![];
    let mut ends = vec![];
    let mut phasesets = vec![];
    let mut exons = vec![];
    let mut bam = if args.input == "-" {
        bam::Reader::from_stdin().expect("\n\nError reading alignments from stdin.\nDid you include the file header with -h?\n\n\n\n")
    } else if args.input.starts_with("s3") || args.input.starts_with("https://") {
        setup_ssl_certificates();
        bam::Reader::from_url(&Url::parse(&args.input).expect("Failed to parse URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::Reader::from_path(&args.input)
            .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n")
    };
    if args.input.ends_with(".cram") & args.reference.is_some() {
        // bam.set_cram_option(htslib::CFR_REQUIRED_FIELDS, htslib::sam_fields_SAM_AUX as i32)
        //     .expect("Failed setting cram options");
        bam.set_reference(
            args.reference
                .as_ref()
                .expect("Failed setting reference for CRAM file"),
        )
        .expect("Failed setting reference for CRAM file");
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

    let min_read_len = args.min_read_len;
    // the match statement below is a bit ugly, but it is the only way to get a closure
    // that closure is used for filtering the reads
    // the closure is different depending on inclusion of unmapped reads (--ubam) and the minimum read length (--min-read-len)
    let filter_closure: Box<dyn Fn(&bam::Record) -> bool> = match (args.ubam, args.min_read_len) {
        (false, 0) => Box::new(|record: &bam::Record| {
            // filter out unmapped, no length filter
            record.flags() & htslib::BAM_FUNMAP as u16 == 0
        }),
        (false, l) if l > 0 => Box::new(|record: &bam::Record| {
            // filter out unmapped, with a length filter
            record.flags() & htslib::BAM_FUNMAP as u16 == 0 && record.seq_len() > min_read_len
        }),
        // keep unmapped reads, no length filter
        (true, 0) => Box::new(|_: &bam::Record| true),
        (true, l) if l > 0 => Box::new(|record: &bam::Record| {
            // only length filter, keep unmapped
            record.seq_len() > min_read_len
        }),
        // the pattern below should be unreachable, as the min_read_len is either zero or positive
        (false, _) | (true, _) => unreachable!(),
    };
    let mut all_counts = 0;
    for read in bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|record| record.flags() & (htslib::BAM_FSECONDARY) as u16 == 0)
        .inspect(|_| all_counts += 1)
        .filter(|read| filter_closure(read))
    {
        let read_length = read.seq_len() as u128 - softclipped_bases(&read);
        lengths.push(read_length);
        if !read.is_supplementary() {
            num_reads += 1;
        }
        if args.karyotype || args.phased {
            tids.push(read.tid());
        }
        if args.phased {
            starts.push(read.pos());
            ends.push(read.reference_end());
            phasesets.push(get_phaseset(&read));
        }
        if args.spliced {
            exons.push(get_exon_number(&read));
        }
        if args.ubam {
            // For unmapped reads, estimate accuracy from per-base Q-scores
            let accuracy = qscore_to_accuracy(&read);
            identities.push(accuracy);
            if hist_requested {
                let phred = crate::utils::accuracy_to_phred(accuracy);
                let index = if phred < 40 { phred } else { 40 };
                q_score_counts[index] += 1;
                q_score_bases[index] += read_length;
            }
        } else {
            let identity = gap_compressed_identity(read);
            identities.push(identity);
            if hist_requested {
                let phred = crate::utils::accuracy_to_phred(identity);
                let index = if phred < 40 { phred } else { 40 };
                q_score_counts[index] += 1;
                q_score_bases[index] += read_length;
            }
        }
    }
    if let Some(s) = &args.arrow {
        match args.ubam {
            true => crate::feather::save_as_arrow_ubam(
                s.to_string(),
                lengths.iter().map(|x| *x as u64).collect(),
                identities.clone(),
            ),
            false => crate::feather::save_as_arrow(
                s.to_string(),
                lengths.iter().map(|x| *x as u64).collect(),
                identities.clone(),
            ),
        }
    }

    // sort vectors in descending order (required for N50/N75)
    lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    identities.par_sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
    (
        Data {
            lengths: Some(lengths),
            num_reads,
            all_counts,
            identities: Some(identities),
            q_score_hist: if hist_requested {
                Some(QScoreHistogramData {
                    counts: q_score_counts,
                    bases: q_score_bases,
                })
            } else {
                None
            },
            tids: if args.karyotype || args.phased {
                Some(tids)
            } else {
                None
            },
            starts: if args.phased { Some(starts) } else { None },
            ends: if args.phased { Some(ends) } else { None },
            phasesets: if args.phased { Some(phasesets) } else { None },
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
            100.0
                * (1.0
                    - ((get_nm_tag(&record) - gap_size + gap_count) as f64
                        / (matches + gap_count) as f64))
        }
    }
}

/// Computes the mean estimated accuracy (%) from per-base Q-scores
/// Q-score is Phred-scaled: Q = -10 * log10(P_error)
/// This function converts each Q-score to probability of correctness
/// and returns the average as a percentage.
fn qscore_to_accuracy(record: &bam::Record) -> f64 {
    let quals = record.qual();
    if quals.is_empty() || quals.iter().all(|&q| q == 255) {
        // 255 indicates missing quality - return 0.0 as fallback
        return 0.0;
    }

    let sum_accuracy: f64 = quals
        .iter()
        .map(|&q| {
            // P_error = 10^(-Q/10), P_correct = 1 - P_error
            1.0 - 10_f64.powf(-(q as f64) / 10.0)
        })
        .sum();

    100.0 * sum_accuracy / quals.len() as f64
}

fn get_nm_tag(record: &bam::Record) -> u32 {
    match record.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => u32::from(v),
            Aux::U16(v) => u32::from(v),
            Aux::U32(v) => v,
            Aux::I8(v) => u32::try_from(v).expect("Identified a negative NM tag"),
            Aux::I16(v) => u32::try_from(v).expect("Identified a negative NM tag"),
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
            Aux::I8(v) => Some(u32::try_from(v).unwrap_or_else(|_| {
                panic!(
                    "Invalid: Identified a negative PS tag at {}",
                    std::str::from_utf8(record.qname())
                        .expect("Failed to convert read name to string")
                )
            })),
            Aux::I16(v) => Some(u32::try_from(v).unwrap_or_else(|_| {
                panic!(
                    "Invalid: Identified a negative PS tag at {}",
                    std::str::from_utf8(record.qname())
                        .expect("Failed to convert read name to string")
                )
            })),
            Aux::I32(v) => Some(u32::try_from(v).unwrap_or_else(|_| {
                panic!(
                    "Invalid: Identified a negative PS tag at {}",
                    std::str::from_utf8(record.qname())
                        .expect("Failed to convert read name to string")
                )
            })),
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

fn softclipped_bases(read: &bam::Record) -> u128 {
    (read.cigar().leading_softclips() + read.cigar().trailing_softclips()) as u128
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

        let accuracy = qscore_to_accuracy(&record);
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

        let accuracy = qscore_to_accuracy(&record);
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

        let accuracy = qscore_to_accuracy(&record);
        // Should return 0.0 for missing quality
        assert!((accuracy - 0.0).abs() < 0.01);
    }
}
