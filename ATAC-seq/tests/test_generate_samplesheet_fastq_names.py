#!/usr/bin/env python3
import csv
import subprocess
import tempfile
from pathlib import Path


PIPELINE_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PIPELINE_DIR / "generate_samplesheet.sh"


def touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


def run_case(names):
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        fq_dir = tmp_path / "fastq"
        for name in names:
            touch(fq_dir / name)
        samplesheet = tmp_path / "samplesheet.csv"
        contrasts = tmp_path / "contrasts.csv"
        subprocess.run(
            [
                "bash",
                str(SCRIPT),
                "--input",
                str(fq_dir),
                "--samplesheet",
                str(samplesheet),
                "--contrasts",
                str(contrasts),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        with samplesheet.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))


def assert_single_pe(names, expected_sample):
    rows = run_case(names)
    assert len(rows) == 1, rows
    row = rows[0]
    assert row["sample"] == expected_sample, row
    assert row["layout"] == "PE", row
    assert "R1" in row["r1"] or "_1" in row["r1"], row
    assert "R2" in row["r2"] or "_2" in row["r2"], row


def test_illumina_s_lane_read_001():
    assert_single_pe(
        ["SampleA_S1_L001_R1_001.fastq.gz", "SampleA_S1_L001_R2_001.fastq.gz"],
        "SampleA",
    )


def test_read_before_lane():
    assert_single_pe(
        ["SampleB_R1_L001.fq.gz", "SampleB_R2_L001.fq.gz"],
        "SampleB",
    )


def test_read_001_without_lane():
    assert_single_pe(
        ["SampleC_R1_001.fastq.gz", "SampleC_R2_001.fastq.gz"],
        "SampleC",
    )


def test_numeric_pair_suffix():
    assert_single_pe(["SampleD_1.fq.gz", "SampleD_2.fq.gz"], "SampleD")


def test_sample_subdirectory_simple_read_names():
    rows = run_case(["SampleE/R1.fastq.gz", "SampleE/R2.fastq.gz"])
    assert len(rows) == 1, rows
    assert rows[0]["sample"] == "SampleE", rows
    assert rows[0]["layout"] == "PE", rows


def test_sample_subdirectory_numeric_read_names():
    rows = run_case(["SampleF/1.fq.gz", "SampleF/2.fq.gz"])
    assert len(rows) == 1, rows
    assert rows[0]["sample"] == "SampleF", rows
    assert rows[0]["layout"] == "PE", rows


if __name__ == "__main__":
    test_illumina_s_lane_read_001()
    test_read_before_lane()
    test_read_001_without_lane()
    test_numeric_pair_suffix()
    test_sample_subdirectory_simple_read_names()
    test_sample_subdirectory_numeric_read_names()
    print("generate_samplesheet FASTQ naming tests OK")
