import subprocess
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import filecmp

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
OUTPUT = ROOT / "tests" / "output"

@pytest.fixture(scope="session", autouse=True)
def setup_data():
    """Download or prepare test data if not already present."""
    files = [
        "IonXpress_020_rawlib.hg19.bam",
        "region.vcf",
        "hg19.fa"
    ]
    if not all((DATA / f).exists() for f in files):
        print("Downloading data...")
        subprocess.run(["bash", str(DATA / "setup.sh")], check=True)

@pytest.fixture(scope="session")
def run_pipeline():
    """Run varlociraptor.sh once and return the process result."""
    print("Running pipeline...")
    cmd = [
        "bash", "varlociraptor.sh",
        "--id", "test",
        "--tumor_cell_quantity", "0.5",
        "--bam", "tests/data/IonXpress_020_rawlib.hg19.bam",
        "--vcf", "tests/data/region.vcf",
        "--reference", "tests/data/hg19.fa",
        "--n_cores", "32",
        "--max_depth", "1000",
        "--output_folder", "tests/output"
    ]
    result = subprocess.run(
        cmd, cwd=ROOT, capture_output=True, text=True
    )
    assert result.returncode == 0, f"Pipeline failed:\n{result.stderr}"
    return result

def compare_files(file1: Path, file2: Path):
    """Helper to assert two files are identical."""
    assert file1.exists(), f"Missing file: {file1}"
    assert file2.exists(), f"Missing file: {file2}"

    # Compare binary content directly (safe for text and binary)
    same = filecmp.cmp(file1, file2, shallow=False)
    if not same:
        raise AssertionError(f"Files differ: {file1} vs {file2}")
    
def compare_tsv(file1: Path, file2: Path, str_cols=None, tol=1e-4):
    if str_cols is None:
        str_cols = []
    assert file1.exists(), f"Missing file: {file1}"
    assert file2.exists(), f"Missing file: {file2}"

    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")

    # Ensure same shape
    assert df1.shape == df2.shape, f"Different shapes: {file1} {df1.shape} vs {file2} {df2.shape}"

    # Identify numeric PROB_ columns
    prob_cols = [c for c in df1.columns if c.startswith("PROB_")]
    cols_to_compare = str_cols + prob_cols

    for col in cols_to_compare:
        assert col in df1.columns, f"Missing column {col} in {file1}"
        assert col in df2.columns, f"Missing column {col} in {file2}"

        if col.startswith("PROB_"):
            # Numeric comparison within tolerance
            diff = np.abs(df1[col] - df2[col])
            assert np.allclose(df1[col], df2[col], atol=tol, equal_nan=True), (
                f"Column {col} differs beyond tolerance {tol} "
                f"in {file1.name} (max diff={diff.max()})"
            )
        else:
            # Exact string match
            assert (df1[col].astype(str) == df2[col].astype(str)).all(), (
                f"String column {col} differs in {file1.name}"
            )

def test_data_files_match(run_pipeline):
    """Check that data/ intermediate results are identical between test and v1 outputs."""
    print("Testing data files match")
    files = [
        "candidates_unnormalized.tsv",
        "candidates_all.tsv",
        "candidates.tsv",
    ]
    for f in files:
        compare_files(OUTPUT / "test" / "data" / f, OUTPUT / "v1" / "data" / f)
    compare_tsv(OUTPUT / "test" / "data" / "calls_filt.tsv", OUTPUT / "v1" / "data" / "calls_filt.tsv")


def test_result_files_match(run_pipeline):
    """Check that result files are identical between test and v1 outputs."""
    compare_tsv(OUTPUT / "test" / "results" / "varlociraptor_probs.tsv", OUTPUT / "v1" / "results" / "varlociraptor_probs.tsv")
