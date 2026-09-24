use std::path::PathBuf;
use std::process::Command;

fn run_cramino_json(args: &[&str]) -> serde_json::Value {
    let output = Command::new(env!("CARGO_BIN_EXE_cramino"))
        .args(args)
        .output()
        .expect("Failed to run cramino");
    assert!(
        output.status.success(),
        "cramino failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("JSON output")
}

fn test_data(name: &str) -> String {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("test-data");
    path.push(name);
    path.to_string_lossy().into_owned()
}

/// The read length statistics of a bam file are the same as those of the ubam with the same reads,
/// while the aligned length statistics only exist meaningfully for the bam file
#[test]
fn bam_and_ubam_have_the_same_read_length_stats() {
    let bam = run_cramino_json(&["--format", "json", &test_data("small-test-phased.bam")]);
    let ubam = run_cramino_json(&[
        "--format",
        "json",
        "--ubam",
        &test_data("small-test-ubam.bam"),
    ]);
    let (bam, ubam) = (&bam["read_stats"], &ubam["read_stats"]);
    for key in ["n50", "n75", "median_length", "mean_length"] {
        assert_eq!(bam[key], ubam[key], "{key} differs between bam and ubam");
    }
    assert_eq!(bam["n50"], 21885);
    assert_eq!(bam["median_length"], 8283.5);
    assert_eq!(bam["n50_aligned"], 19294);
    assert_eq!(bam["median_length_aligned"], 6767.5);
}

#[test]
fn ubam_has_no_mean_coverage() {
    let ubam = run_cramino_json(&[
        "--format",
        "json",
        "--ubam",
        &test_data("small-test-ubam.bam"),
    ]);
    assert!(ubam["read_stats"]["mean_coverage"].is_null());
}

#[test]
fn no_reads_passing_filters_have_zero_coverage() {
    // with reference sequences in the header, the coverage is 0 rather than unavailable
    let bam = run_cramino_json(&[
        "--format",
        "json",
        "--min-read-len",
        "10000000",
        &test_data("small-test-phased.bam"),
    ]);
    assert_eq!(bam["read_stats"]["mean_coverage"], 0.0);
}
