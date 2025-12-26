use std::path::PathBuf;
use std::process::Command;

fn run_cramino_json(args: Vec<String>) -> serde_json::Value {
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

fn test_bam_path() -> String {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("test-data");
    path.push("small-test-phased.bam");
    path.to_string_lossy().into_owned()
}

fn test_hist_path() -> String {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("json-hist.txt");
    path.to_string_lossy().into_owned()
}

#[test]
fn json_omits_histograms_without_hist_flags() {
    let bam_path = test_bam_path();
    let args = vec!["--format".to_string(), "json".to_string(), bam_path];
    let json_value = run_cramino_json(args);
    assert!(json_value.get("histograms").is_none());
}

#[test]
fn json_includes_histograms_with_hist_flag() {
    let bam_path = test_bam_path();
    let hist_arg = format!("--hist={}", test_hist_path());
    let args = vec![
        "--format".to_string(),
        "json".to_string(),
        hist_arg,
        bam_path,
    ];
    let json_value = run_cramino_json(args);
    assert!(json_value.get("histograms").is_some());
    assert!(json_value["histograms"]["read_length"]["bins"].is_array());
    assert!(json_value["histograms"]["q_score"]["bins"].is_array());
}
