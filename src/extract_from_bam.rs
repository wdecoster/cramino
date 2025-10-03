use bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};
use std::env;
use url::Url;
use rayon::prelude::*;
use log::warn;


pub struct Data {
    pub lengths: Option<Vec<u128>>,
    pub num_reads: usize,
    pub all_counts: usize,
    pub identities: Option<Vec<f64>>,
    pub tids: Option<Vec<i32>>,
    pub starts: Option<Vec<i64>>,
    pub ends: Option<Vec<i64>>,
    pub phasesets: Option<Vec<Option<u32>>>,
    pub exons: Option<Vec<usize>>,
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
            env::set_var("CURL_CA_BUNDLE", path);
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
            record.flags() & htslib::BAM_FUNMAP as u16 == 0
                && record.seq_len() > min_read_len
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
        lengths.push(read.seq_len() as u128 - softclipped_bases(&read));
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
        if !args.ubam {
            identities.push(gap_compressed_identity(read));
        }
    }
    if let Some(s) = &args.arrow {
        match args.ubam {
            true => crate::feather::save_as_arrow_ubam(
                s.to_string(),
                lengths.iter().map(|x| *x as u64).collect(),
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
            identities: if !args.ubam { Some(identities) } else { None },
            tids: if args.karyotype || args.phased {
                Some(tids)
            } else {
                None
            },
            starts: if args.phased { Some(starts) } else { None },
            ends: if args.phased { Some(ends) } else { None },
            phasesets: if args.phased { Some(phasesets) } else { None },
            exons: if args.spliced { Some(exons) } else { None },
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
