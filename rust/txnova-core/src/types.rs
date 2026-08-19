use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignerFamily {
    Star,
    Hisat2,
    Bowtie2,
    Minimap2,
    Other,
}

impl AlignerFamily {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Star => "STAR",
            Self::Hisat2 => "HISAT2",
            Self::Bowtie2 => "Bowtie2",
            Self::Minimap2 => "minimap2",
            Self::Other => "other",
        }
    }

    pub fn accepted(self) -> bool {
        matches!(self, Self::Star | Self::Hisat2)
    }

    pub fn forbidden(self) -> bool {
        matches!(self, Self::Bowtie2 | Self::Minimap2)
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct SampleIn {
    pub sample_id: String,
    pub bam: String,
    #[serde(default)]
    pub group: String,
    pub strandedness: String,
    pub replicate: u32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PreflightInput {
    pub assembly: String,
    pub library_layout: String,
    pub de_enabled: bool,
    pub treat_min_detected_replicates: u32,
    #[serde(default)]
    pub threads: usize,
    pub samples: Vec<SampleIn>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SampleOut {
    pub sample_id: String,
    pub n_sq: usize,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct PreflightReport {
    pub ok: bool,
    pub assembly: String,
    pub sq_sha256: String,
    pub aligner_family: String,
    pub strandedness: String,
    pub library_layout: String,
    pub n_control: usize,
    pub n_treat: usize,
    pub dup_flag_seen: serde_json::Map<String, serde_json::Value>,
    pub gtf_seqnames: usize,
    pub samples: Vec<SampleOut>,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
    pub threads: ThreadPlan,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct ThreadPlan {
    pub requested: usize,
    pub cpus: usize,
    pub usable: usize,
    pub available_bytes: Option<u64>,
    pub bam_workers: usize,
    pub auto: bool,
}
