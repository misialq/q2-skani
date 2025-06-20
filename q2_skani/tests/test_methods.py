import tempfile
from pathlib import Path
from typing import List
from unittest.mock import patch, MagicMock
import subprocess

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_skani.skani import (
    _construct_triangle_cmd,
    _run_skani,
    _process_skani_matrix,
    compare_seqs,
)


class SkaniMethodTests(TestPluginBase):
    package = "q2_skani.tests"

    def setUp(self):
        super().setUp()
        self.test_matrix = self.get_data_path("test.matrix")
        self.invalid_matrix = self.get_data_path("test_invalid.matrix")
        self.na_matrix = self.get_data_path("test_na.matrix")

    # ------------------------------------------------------------------
    # _construct_triangle_cmd
    # ------------------------------------------------------------------
    def test_construct_triangle_cmd_minimal(self):
        cmd = _construct_triangle_cmd(
            fasta_list="genomes.txt",
            output_file="out.tsv",
            skani_args={},
        )
        exp = [
            "skani",
            "triangle",
            "-v",
            "--distance",
            "-l",
            "genomes.txt",
            "-o",
            "out.tsv",
        ]
        self.assertEqual(cmd, exp)

    def test_construct_triangle_cmd_all_opts(self):
        cmd = _construct_triangle_cmd(
            fasta_list="genomes.txt",
            output_file="out.tsv",
            skani_args={
                "threads": 8,
                "min_af": 10.0,
                "compression": 50,
                "marker_c": 500,
                "screen": 90.0,
                "ci": True,
                "detailed": True,
                "diagonal": True,
                "sparse": True,
                "full_matrix": True,
                "median": True,
                "no_learned_ani": True,
                "robust": True,
                "faster_small": True,
                "preset": "fast",
            },
        )
        exp = [
            "skani",
            "triangle",
            "-v",
            "--distance",
            "-l",
            "genomes.txt",
            "-o",
            "out.tsv",
            "-t",
            "8",
            "--min-af",
            "10.0",
            "-c",
            "50",
            "-m",
            "500",
            "-s",
            "90.0",
            "--ci",
            "--detailed",
            "--diagonal",
            "--sparse",
            "--full-matrix",
            "--median",
            "--no-learned-ani",
            "--robust",
            "--faster-small",
            "--fast",
        ]
        self.assertEqual(cmd, exp)

    # ------------------------------------------------------------------
    # _run_skani
    # ------------------------------------------------------------------
    def test_run_skani_success(self):
        with patch("subprocess.run") as run:
            run.return_value = MagicMock()
            _run_skani(["skani"])
            run.assert_called_once_with(
                ["skani"],
                check=True,
                capture_output=True,
                text=True,
            )

    def test_run_skani_failure(self):
        error = subprocess.CalledProcessError(1, ["skani"], stdout="out", stderr="err")
        with patch("subprocess.run", side_effect=error):
            with self.assertRaises(RuntimeError) as cm:
                _run_skani(["skani"])
        self.assertIn("Skani failed with exit code 1", str(cm.exception))
        self.assertIn("stdout:\nout", str(cm.exception))
        self.assertIn("stderr:\nerr", str(cm.exception))

    # ------------------------------------------------------------------
    # _process_skani_matrix
    # ------------------------------------------------------------------
    def test_process_skani_matrix(self):
        df = _process_skani_matrix(self.test_matrix)
        self.assertEqual(df.shape, (3, 3))
        self.assertEqual(list(df.index), ["genome1", "genome2", "genome3"])
        self.assertEqual(df.loc["genome1", "genome2"], 0.1)

    def test_process_skani_matrix_invalid(self):
        with self.assertRaises(Exception):
            _process_skani_matrix(self.invalid_matrix)

    # ------------------------------------------------------------------
    # compare_seqs
    # ------------------------------------------------------------------
    def test_compare_seqs(self):
        mags = MockMAGSequencesDirFmt(["g1.fasta", "g2.fasta"])

        def fake_run(cmd):
            out_file = cmd[cmd.index("-o") + 1]
            Path(out_file).write_text(Path(self.test_matrix).read_text())

        with patch("q2_skani.skani._run_skani", side_effect=fake_run):
            dm = compare_seqs(mags)

        self.assertEqual(dm.shape, (3, 3))
        self.assertEqual(dm.ids, ["genome1", "genome2", "genome3"])


class MockMAGSequencesDirFmt:
    def __init__(self, fasta_files: List[str]):
        self.path = Path(tempfile.mkdtemp())
        for f in fasta_files:
            (self.path / f).touch()
