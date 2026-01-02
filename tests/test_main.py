"""
Tests for main function via subprocess.
"""

from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from ht1.ht1 import ErrorCode

TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_FILES = list(TEST_DATA_DIR.glob("**/SRA/*.raw"))


def run_main(*args):
    """Run search.py as subprocess with given arguments."""
    result = subprocess.run(
        [sys.executable, "-m", "ht1.ht1", *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return result


# --- Valid file tests ---


@pytest.mark.parametrize("filepath", TEST_FILES, ids=lambda p: p.parent.parent.name)
def test_main_completes_on_valid_files(filepath, tmp_path):
    """Main completes without error on all valid test files."""
    # Copy file to tmp_path
    tmp_file = tmp_path / filepath.name
    shutil.copy(filepath, tmp_file)

    result = run_main(str(tmp_file))

    assert result.returncode == ErrorCode.OK, f"stderr: {result.stderr}"


@pytest.mark.parametrize("filepath", TEST_FILES[:3], ids=lambda p: p.parent.parent.name)
def test_main_output_location(filepath, tmp_path):
    """Output file is created in same directory as input with correct name."""
    tmp_file = tmp_path / filepath.name
    shutil.copy(filepath, tmp_file)

    run_main(str(tmp_file), "--threshold", "0.1")  # low threshold to ensure hits

    expected_output = tmp_path / f"{tmp_file.stem}_trigger.txt"
    assert expected_output.exists()


def test_main_output_format(tmp_path):
    """Output file has expected format: nhits start end."""
    # Use first available test file
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    run_main(str(tmp_file), "--threshold", "0.1")

    output_file = tmp_path / f"{tmp_file.stem}_trigger.txt"
    if output_file.exists():
        content = output_file.read_text()
        parts = content.split()
        assert len(parts) == 3
        # All parts should be integers
        nhits, start, end = int(parts[0]), int(parts[1]), int(parts[2])
        assert nhits > 0
        assert start < end


def test_main_no_output_when_no_hits(tmp_path):
    """No output file created when there are no hits (high threshold)."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    # Very high threshold should produce no hits
    run_main(str(tmp_file), "--threshold", "100.0")

    output_file = tmp_path / f"{tmp_file.stem}_trigger.txt"
    assert not output_file.exists()


# --- Invalid file tests ---


def test_main_file_not_found(tmp_path):
    """Non-existent file returns INVALID_FILE error code."""
    result = run_main(str(tmp_path / "nonexistent.raw"))

    assert result.returncode == ErrorCode.INVALID_FILE
    assert "does not exist" in result.stdout


def test_main_invalid_sra(tmp_path):
    """Invalid SRA file returns INVALID_FILE error code."""
    bad_file = tmp_path / "bad.raw"
    bad_file.write_bytes(b"too short")

    result = run_main(str(bad_file))

    assert result.returncode == ErrorCode.INVALID_FILE
    assert "not valid" in result.stdout


# --- Invalid parameter tests ---


def test_main_invalid_size_zero(tmp_path):
    """--size 0 returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--size", "0")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS
    assert "--size" in result.stdout


def test_main_invalid_size_negative(tmp_path):
    """--size -1 returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--size", "-1")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS


def test_main_invalid_maxtest_not_power_of_2(tmp_path):
    """--maxtest 3 (not power of 2) returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--maxtest", "3")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS
    assert "--maxtest" in result.stdout


def test_main_invalid_maxtest_zero(tmp_path):
    """--maxtest 0 returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--maxtest", "0")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS


def test_main_invalid_threshold_zero(tmp_path):
    """--threshold 0 returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--threshold", "0")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS
    assert "--threshold" in result.stdout


def test_main_invalid_threshold_negative(tmp_path):
    """--threshold -5 returns INVALID_PARAMETERS error code."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--threshold", "-5")

    assert result.returncode == ErrorCode.INVALID_PARAMETERS


# --- Custom parameter tests ---


def test_main_custom_parameters(tmp_path):
    """Main accepts custom size, maxtest, and threshold."""
    if not TEST_FILES:
        pytest.skip("No test files available")

    tmp_file = tmp_path / TEST_FILES[0].name
    shutil.copy(TEST_FILES[0], tmp_file)

    result = run_main(str(tmp_file), "--size", "100", "--maxtest", "16", "--threshold", "4.0")

    assert result.returncode == ErrorCode.OK


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
