use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Metrics {
    pub file_info: FileInfo,
    pub alignment_stats: AlignmentStats,
    pub read_stats: ReadStats,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub identity_stats: Option<IdentityStats>,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_stats: Option<PhaseStats>,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub karyotype_stats: Option<Vec<ChromosomeData>>,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub splice_stats: Option<SpliceStats>,
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

#[derive(Serialize, Deserialize, Debug)]
pub struct ReadStats {
    pub yield_gb: f64,
    pub mean_coverage: f64,
    pub yield_gb_long: f64,
    pub n50: u128,
    pub n75: u128,
    pub median_length: f64,
    pub mean_length: f64,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct IdentityStats {
    pub median_identity: f64,
    pub mean_identity: f64,
    pub modal_identity: f64,
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
    pub count: usize,
    pub normalized_count: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct SpliceStats {
    pub median_exons: usize,
    pub mean_exons: f32,
    pub fraction_unspliced: f32,
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
            read_stats: ReadStats {
                yield_gb: 0.0,
                mean_coverage: 0.0,
                yield_gb_long: 0.0,
                n50: 0,
                n75: 0,
                median_length: 0.0,
                mean_length: 0.0,
            },
            identity_stats: None,
            phase_stats: None,
            karyotype_stats: None,
            splice_stats: None,
        }
    }
}