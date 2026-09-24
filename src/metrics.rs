use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Metrics {
    pub file_info: FileInfo,
    pub alignment_stats: AlignmentStats,
    pub read_stats: ReadStats,

    /// Gap-compressed identity of aligned reads
    #[serde(skip_serializing_if = "Option::is_none")]
    pub identity_stats: Option<IdentityStats>,

    /// Identity estimated from the base qualities, with --ubam
    #[serde(skip_serializing_if = "Option::is_none")]
    pub estimated_identity_stats: Option<EstimatedIdentityStats>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_stats: Option<PhaseStats>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub karyotype_stats: Option<Vec<ChromosomeData>>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub splice_stats: Option<SpliceStats>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub histograms: Option<Histograms>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FileInfo {
    pub name: String,
    pub path: String,
    pub creation_time: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct AlignmentStats {
    pub num_alignments: usize,
    pub percent_from_total: f64,
    pub num_reads: usize,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct ReadStats {
    /// Yield and coverage are of the aligned bases
    pub yield_gb: f64,
    /// Not available without reference sequences in the header, e.g. for a ubam
    pub mean_coverage: Option<f64>,
    pub yield_gb_long: f64,
    /// Statistics of the full read lengths, i.e. of the primary alignments including clipped bases
    pub n50: u128,
    pub n75: u128,
    pub median_length: f64,
    pub mean_length: f64,
    /// Statistics of the aligned lengths of the primary and supplementary alignments, without clipped bases
    pub n50_aligned: u128,
    pub n75_aligned: u128,
    pub median_length_aligned: f64,
    pub mean_length_aligned: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct IdentityStats {
    pub median_identity: f64,
    pub mean_identity: f64,
    pub modal_identity: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct EstimatedIdentityStats {
    pub median_estimated_identity: f64,
    pub mean_estimated_identity: f64,
    pub modal_estimated_identity: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PhaseStats {
    pub fraction_phased: f32,
    pub num_phaseblocks: usize,
    pub total_bases_phased_gb: f64,
    pub median_phaseblock_length: f64,
    pub n50_phaseblock_length: i64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ChromosomeData {
    pub chromosome: String,
    /// Number of reads, i.e. primary alignments
    pub count: usize,
    /// Reads per bp, relative to the median over all chromosomes with reads
    pub normalized_count: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SpliceStats {
    pub median_exons: usize,
    pub mean_exons: f32,
    pub fraction_unspliced: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Histograms {
    pub read_length: Histogram,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub q_score: Option<Histogram>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Histogram {
    pub step: u64,
    pub max_value: u64,
    pub bins: Vec<HistogramBin>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct HistogramBin {
    pub start: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<u64>,
    pub count: u64,
    pub bases: u128,
}

impl Metrics {
    pub fn new(file_info: FileInfo) -> Self {
        Metrics {
            file_info,
            alignment_stats: AlignmentStats {
                num_alignments: 0,
                percent_from_total: 0.0,
                num_reads: 0,
            },
            read_stats: ReadStats::default(),
            identity_stats: None,
            estimated_identity_stats: None,
            phase_stats: None,
            karyotype_stats: None,
            splice_stats: None,
            histograms: None,
        }
    }
}
