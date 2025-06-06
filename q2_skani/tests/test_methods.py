# ----------------------------------------------------------------------------
# Copyright (c) 2025, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
from pathlib import Path
from typing import List
from unittest.mock import patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_skani.skani import (
    compare_seqs,
    _construct_skani_cmd,
    _process_skani_matrix,
    compare_seqs_skani,
)


class FastANITests(TestPluginBase):
    package = "q2_skani.tests"

    def setUp(self):
        """Set up test data."""
        super().setUp()
        self.test_data_dir = Path(self.temp_dir.name) / "test_data"
        self.test_data_dir.mkdir()
        # Create a mock for subprocess.run that will be used in all tests
        self.mock_subprocess = patch("subprocess.run").start()
        self.addCleanup(patch.stopall)

    def test_fastani_never_called(self):
        """Test that FastANI is never actually called."""
        # Create mock MAG sequences
        mags = MockMAGSequencesDirFmt(["genome1.fasta", "genome2.fasta"])

        # Mock the FastANI command execution
        def mock_run_fastani(cmd: List[str]):
            # Create mock output files
            output_dir = Path(cmd[cmd.index("--output") + 1]).parent
            matrix_file = output_dir / "fastani_output.tsv.matrix"
            # Copy the test matrix file to the output location
            with open(self.get_data_path("test.matrix"), "r") as src:
                with open(matrix_file, "w") as dst:
                    dst.write(src.read())

        self.monkeypatch.setattr("q2_skani._methods._run_fastani", mock_run_fastani)

        # Run the comparison
        compare_seqs(
            genomes=mags,
            kmer=16,
            fragLen=3000,
            threads=1,
            minFraction=0.2,
            maxRatioDiff=0.05,
        )

        # Verify that subprocess.run was never called
        self.mock_subprocess.assert_not_called()

    def test_construct_fastani_cmd(self):
        """Test command construction with various parameters."""
        # Test with minimal required parameters
        cmd = _construct_fastani_cmd(
            query_list="query.txt",
            ref_list="ref.txt",
            output_file="output.tsv",
            fastani_args={},
        )
        self.assertEqual(
            cmd,
            [
                "fastANI",
                "--queryList",
                "query.txt",
                "--refList",
                "ref.txt",
                "--output",
                "output.tsv",
                "--matrix",
            ],
        )

        # Test with all optional parameters
        cmd = _construct_fastani_cmd(
            query_list="query.txt",
            ref_list="ref.txt",
            output_file="output.tsv",
            fastani_args={
                "kmer": 16,
                "fragLen": 3000,
                "threads": 4,
                "minFraction": 0.2,
                "maxRatioDiff": 0.05,
            },
        )
        self.assertEqual(
            cmd,
            [
                "fastANI",
                "--queryList",
                "query.txt",
                "--refList",
                "ref.txt",
                "--output",
                "output.tsv",
                "--matrix",
                "--kmer",
                "16",
                "--fragLen",
                "3000",
                "--threads",
                "4",
                "--minFraction",
                "0.2",
                "--maxRatioDiff",
                "0.05",
            ],
        )

    def test_process_fastani_matrix(self):
        """Test processing of FastANI matrix output."""
        # Use the test matrix file
        matrix_file = self.get_data_path("test.matrix")

        # Process the matrix
        df = _process_fastani_matrix(str(matrix_file))

        # Check the results
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape, (3, 3))
        self.assertEqual(list(df.index), ["genome1", "genome2", "genome3"])
        self.assertEqual(list(df.columns), ["genome1", "genome2", "genome3"])
        self.assertEqual(df.loc["genome1", "genome2"], 99.9)  # 100 - 0.1
        self.assertEqual(df.loc["genome2", "genome1"], 99.9)  # Symmetric
        self.assertEqual(df.loc["genome1", "genome1"], 0.0)  # Self-comparison

    def test_process_fastani_matrix_with_na(self):
        """Test processing of FastANI matrix with NA values."""
        # Use the test matrix file with NA values
        matrix_file = self.get_data_path("test_na.matrix")

        # Process the matrix
        df = _process_fastani_matrix(str(matrix_file))

        # Check the results
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape, (3, 3))
        self.assertTrue(pd.isna(df.loc["genome1", "genome2"]))
        self.assertTrue(pd.isna(df.loc["genome2", "genome1"]))

    def test_process_fastani_matrix_invalid_file(self):
        """Test error handling for invalid matrix file."""
        # Use the invalid test matrix file
        matrix_file = self.get_data_path("test_invalid.matrix")

        # Test that processing raises an error
        with self.assertRaises(ValueError):
            _process_fastani_matrix(str(matrix_file))

    def test_compare_seqs(self):
        """Test the main compare_seqs function."""
        # Create mock MAG sequences
        mags = MockMAGSequencesDirFmt(["genome1.fasta", "genome2.fasta"])

        # Mock the FastANI command execution
        def mock_run_fastani(cmd: List[str]):
            # Create mock output files
            output_dir = Path(cmd[cmd.index("--output") + 1]).parent
            matrix_file = output_dir / "fastani_output.tsv.matrix"
            # Copy the test matrix file to the output location
            with open(self.get_data_path("test.matrix"), "r") as src:
                with open(matrix_file, "w") as dst:
                    dst.write(src.read())

        self.monkeypatch.setattr("q2_skani._methods._run_fastani", mock_run_fastani)

        # Run the comparison
        result = compare_seqs(
            genomes=mags,
            kmer=16,
            fragLen=3000,
            threads=1,
            minFraction=0.2,
            maxRatioDiff=0.05,
        )

        # Check the results
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (3, 3))
        self.assertEqual(list(result.index), ["genome1", "genome2", "genome3"])
        self.assertEqual(result.loc["genome1", "genome2"], 99.9)

    def test_compare_seqs_no_fasta_files(self):
        """Test error handling when no FASTA files are found."""
        # Create mock MAG sequences with no FASTA files
        mags = MockMAGSequencesDirFmt([])

        # Test that comparison raises an error
        with self.assertRaisesRegex(RuntimeError, "No FASTA files found"):
            compare_seqs(genomes=mags)

    def test_compare_seqs_fastani_error(self):
        """Test error handling when FastANI fails."""
        # Create mock MAG sequences
        mags = MockMAGSequencesDirFmt(["genome1.fasta"])

        # Mock FastANI to raise an error
        def mock_run_fastani(cmd: List[str]):
            raise RuntimeError("FastANI failed")

        self.monkeypatch.setattr("q2_skani._methods._run_fastani", mock_run_fastani)

        # Test that comparison raises an error
        with self.assertRaisesRegex(RuntimeError, "Failed to run FastANI comparison"):
            compare_seqs(genomes=mags)

    def test_construct_skani_cmd(self):
        """Test skani command construction with various parameters."""
        # Test with minimal required parameters
        cmd = _construct_skani_cmd(
            query_list="query.txt",
            ref_list="ref.txt",
            output_file="output.tsv",
            skani_args={},
        )
        self.assertEqual(
            cmd,
            [
                "skani",
                "dist",
                "-q",
                "query.txt",
                "-r",
                "ref.txt",
                "-o",
                "output.tsv",
            ],
        )

        # Test with all optional parameters
        cmd = _construct_skani_cmd(
            query_list="query.txt",
            ref_list="ref.txt",
            output_file="output.tsv",
            skani_args={
                "threads": 4,
                "min_fraction": 0.2,
                "min_ani": 80.0,
                "small_genomes": True,
            },
        )
        self.assertEqual(
            cmd,
            [
                "skani",
                "dist",
                "-q",
                "query.txt",
                "-r",
                "ref.txt",
                "-o",
                "output.tsv",
                "-t",
                "4",
                "-m",
                "0.2",
                "-c",
                "80.0",
                "--small-genomes",
            ],
        )

    def test_process_skani_matrix(self):
        """Test processing of skani matrix output."""
        # Use the test matrix file
        matrix_file = self.get_data_path("test.matrix")

        # Process the matrix
        df = _process_skani_matrix(str(matrix_file))

        # Check the results
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape, (3, 3))
        self.assertEqual(list(df.index), ["genome1", "genome2", "genome3"])
        self.assertEqual(list(df.columns), ["genome1", "genome2", "genome3"])
        self.assertEqual(df.loc["genome1", "genome2"], 99.9)  # 100 - 0.1
        self.assertEqual(df.loc["genome2", "genome1"], 99.9)  # Symmetric
        self.assertEqual(df.loc["genome1", "genome1"], 0.0)  # Self-comparison

    def test_process_skani_matrix_with_na(self):
        """Test processing of skani matrix with NA values."""
        # Use the test matrix file with NA values
        matrix_file = self.get_data_path("test_na.matrix")

        # Process the matrix
        df = _process_skani_matrix(str(matrix_file))

        # Check the results
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape, (3, 3))
        self.assertTrue(pd.isna(df.loc["genome1", "genome2"]))
        self.assertTrue(pd.isna(df.loc["genome2", "genome1"]))

    def test_process_skani_matrix_invalid_file(self):
        """Test error handling for invalid matrix file."""
        # Use the invalid test matrix file
        matrix_file = self.get_data_path("test_invalid.matrix")

        # Test that processing raises an error
        with self.assertRaises(ValueError):
            _process_skani_matrix(str(matrix_file))

    def test_compare_seqs_skani(self):
        """Test the main compare_seqs_skani function."""
        # Create mock MAG sequences
        mags = MockMAGSequencesDirFmt(["genome1.fasta", "genome2.fasta"])

        # Mock the skani command execution
        def mock_run_skani(cmd: List[str]):
            # Create mock output files
            output_dir = Path(cmd[cmd.index("-o") + 1]).parent
            matrix_file = output_dir / "skani_output.tsv.matrix"
            # Copy the test matrix file to the output location
            with open(self.get_data_path("test.matrix"), "r") as src:
                with open(matrix_file, "w") as dst:
                    dst.write(src.read())

        self.monkeypatch.setattr("q2_skani._methods._run_skani", mock_run_skani)

        # Run the comparison
        result = compare_seqs_skani(
            genomes=mags,
            threads=4,
            min_fraction=0.2,
            min_ani=80.0,
            small_genomes=True,
        )

        # Check the results
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (3, 3))
        self.assertEqual(list(result.index), ["genome1", "genome2", "genome3"])
        self.assertEqual(result.loc["genome1", "genome2"], 99.9)

    def test_compare_seqs_skani_no_fasta_files(self):
        """Test error handling when no FASTA files are found."""
        # Create mock MAG sequences with no FASTA files
        mags = MockMAGSequencesDirFmt([])

        # Test that comparison raises an error
        with self.assertRaisesRegex(RuntimeError, "No FASTA files found"):
            compare_seqs_skani(genomes=mags)

    def test_compare_seqs_skani_error(self):
        """Test error handling when skani fails."""
        # Create mock MAG sequences
        mags = MockMAGSequencesDirFmt(["genome1.fasta"])

        # Mock skani to raise an error
        def mock_run_skani(cmd: List[str]):
            raise RuntimeError("skani failed")

        self.monkeypatch.setattr("q2_skani._methods._run_skani", mock_run_skani)

        # Test that comparison raises an error
        with self.assertRaisesRegex(RuntimeError, "Failed to run skani comparison"):
            compare_seqs_skani(genomes=mags)


class MockMAGSequencesDirFmt:
    """Mock class for MAGSequencesDirFmt."""

    def __init__(self, fasta_files: List[str]):
        self.path = Path(tempfile.mkdtemp())
        for file in fasta_files:
            (self.path / file).touch()
