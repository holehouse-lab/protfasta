"""
Tests for the pfasta command-line interface (protfasta.scripts.pfasta).

Covers:
- Helper functions (stdout, exit_error, validate, validate_int,
  print_statistical_summary)
- CLI argument parsing and end-to-end invocation of main()
"""

from __future__ import annotations

import os
import sys

import pytest

import protfasta
from protfasta.scripts.pfasta import (
    exit_error,
    main,
    print_statistical_summary,
    stdout,
    validate,
    validate_int,
)
from protfasta import __version__ as VERSION_MAJ

# ── Paths to shared test data ──────────────────────────────────────────────

TEST_DATA_DIR = protfasta._get_data("test_data")
SIMPLE_FILE = os.path.join(TEST_DATA_DIR, "testset_1.fasta")
BADCHAR_FILE = os.path.join(TEST_DATA_DIR, "test_data_with_bad_chars.fa")
NONSTANDARD_FILE = os.path.join(TEST_DATA_DIR, "test_data_with_nonstandard_chars.fa")
DUPLICATE_RECORD_FILE = os.path.join(TEST_DATA_DIR, "testset_duplicate.fasta")
DUPLICATE_SEQ_FILE = os.path.join(TEST_DATA_DIR, "testset_duplicate_seqs.fasta")


def _run_main(*cli_args: str, monkeypatch) -> None:
    """Invoke main() after patching sys.argv with the supplied CLI args."""
    monkeypatch.setattr(sys, "argv", ["pfasta", *cli_args])
    main()


# ── Helper function tests ──────────────────────────────────────────────────


class TestStdout:
    """Tests for the stdout() helper."""

    def test_prints_when_not_silent(self, capsys):
        stdout("hello", False)
        assert capsys.readouterr().out == "hello\n"

    def test_silent_suppresses_output(self, capsys):
        stdout("hello", True)
        assert capsys.readouterr().out == ""


class TestExitError:
    """Tests for the exit_error() helper."""

    def test_exit_error_raises_system_exit(self):
        with pytest.raises(SystemExit):
            exit_error("boom")

    def test_exit_error_prints_message(self, capsys):
        with pytest.raises(SystemExit):
            exit_error("something broke")
        assert "something broke" in capsys.readouterr().out


class TestValidate:
    """Tests for the validate() helper."""

    def test_valid_option_returned_lowercase(self):
        assert validate("FAIL", ["ignore", "fail", "remove"]) == "fail"

    def test_invalid_option_exits(self):
        with pytest.raises(SystemExit):
            validate("nope", ["ignore", "fail", "remove"])

    def test_mixed_case(self):
        assert validate("Remove", ["ignore", "fail", "remove"]) == "remove"


class TestValidateInt:
    """Tests for the validate_int() helper."""

    def test_valid_int(self):
        assert validate_int("10", 0, "--test") == 10

    def test_value_at_minimum(self):
        assert validate_int("5", 5, "--test") == 5

    def test_below_minimum_exits(self):
        with pytest.raises(SystemExit):
            validate_int("3", 5, "--test")

    def test_non_numeric_exits(self):
        with pytest.raises(SystemExit):
            validate_int("abc", 0, "--test")


class TestPrintStatisticalSummary:
    """Tests for the print_statistical_summary() helper."""

    def test_prints_stats(self, capsys):
        data = [["h1", "ACDE"], ["h2", "ACDEFGHIKLMNPQ"], ["h3", "ACD"]]
        print_statistical_summary(data)
        out = capsys.readouterr().out
        assert "[STATS]" in out
        assert "Total number of sequences" in out


# ── CLI integration tests via in-process main() ───────────────────────────


class TestCLIVersion:
    """--version flag."""

    def test_version_flag(self, monkeypatch, capsys):
        with pytest.raises(SystemExit) as exc:
            _run_main("--version", monkeypatch=monkeypatch)
        assert exc.value.code == 0
        assert VERSION_MAJ in capsys.readouterr().out


class TestCLINoArgs:
    """Calling with no arguments should error."""

    def test_no_args_fails(self, monkeypatch):
        with pytest.raises(SystemExit) as exc:
            _run_main(monkeypatch=monkeypatch)
        assert exc.value.code != 0


class TestCLIMissingFile:
    """Pointing at a non-existent file should error."""

    def test_missing_file(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main("no_such_file.fasta", "--silent", monkeypatch=monkeypatch)


class TestCLIBasicRead:
    """Basic invocation with a valid file."""

    def test_basic_read_creates_output(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(SIMPLE_FILE, "-o", outfile, "--silent", monkeypatch=monkeypatch)
        assert os.path.exists(outfile)
        seqs = protfasta.read_fasta(outfile)
        assert len(seqs) > 0

    def test_silent_suppresses_stdout(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        _run_main(SIMPLE_FILE, "-o", outfile, "--silent", monkeypatch=monkeypatch)
        assert capsys.readouterr().out == ""

    def test_verbose_produces_stdout(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        _run_main(SIMPLE_FILE, "-o", outfile, monkeypatch=monkeypatch)
        assert "pfasta version" in capsys.readouterr().out


class TestCLINoOutputFile:
    """--no-outputfile prevents file creation."""

    def test_no_output_file_flag(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _run_main(SIMPLE_FILE, "--no-outputfile", "--silent", monkeypatch=monkeypatch)
        assert not os.path.exists(tmp_path / "output.fasta")


class TestCLIDuplicateRecord:
    """--duplicate-record option."""

    def test_duplicate_record_remove(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            DUPLICATE_RECORD_FILE, "-o", outfile,
            "--duplicate-record", "remove", "--non-unique-header", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)

    def test_duplicate_record_fail(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(Exception):
            _run_main(
                DUPLICATE_RECORD_FILE, "-o", outfile,
                "--duplicate-record", "fail", "--silent",
                monkeypatch=monkeypatch,
            )

    def test_duplicate_record_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            DUPLICATE_RECORD_FILE, "-o", outfile,
            "--duplicate-record", "ignore", "--non-unique-header", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)


class TestCLIDuplicateSequence:
    """--duplicate-sequence option."""

    def test_duplicate_sequence_remove(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            DUPLICATE_SEQ_FILE, "-o", outfile,
            "--duplicate-sequence", "remove", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)

    def test_duplicate_sequence_fail(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(Exception):
            _run_main(
                DUPLICATE_SEQ_FILE, "-o", outfile,
                "--duplicate-sequence", "fail", "--silent",
                monkeypatch=monkeypatch,
            )


class TestCLIInvalidSequence:
    """--invalid-sequence option."""

    def test_invalid_sequence_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            NONSTANDARD_FILE, "-o", outfile,
            "--invalid-sequence", "ignore", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)

    def test_invalid_sequence_fail(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(Exception):
            _run_main(
                NONSTANDARD_FILE, "-o", outfile,
                "--invalid-sequence", "fail", "--silent",
                monkeypatch=monkeypatch,
            )

    def test_invalid_sequence_remove(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        # All sequences in this file have non-standard chars, so removing
        # invalid sequences leaves 0 entries and the CLI exits cleanly.
        with pytest.raises(SystemExit) as exc:
            _run_main(
                NONSTANDARD_FILE, "-o", outfile,
                "--invalid-sequence", "remove", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code == 0

    def test_invalid_sequence_convert_all_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            NONSTANDARD_FILE, "-o", outfile,
            "--invalid-sequence", "convert-all-ignore", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)

    def test_invalid_sequence_convert_res_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            NONSTANDARD_FILE, "-o", outfile,
            "--invalid-sequence", "convert-res-ignore", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)


class TestCLILengthFilters:
    """--shortest-seq and --longest-seq options."""

    def test_shortest_seq(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--shortest-seq", "500", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)
        seqs = protfasta.read_fasta(outfile)
        for seq in seqs.values():
            assert len(seq) > 500

    def test_longest_seq(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--longest-seq", "400", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)
        seqs = protfasta.read_fasta(outfile)
        for seq in seqs.values():
            assert len(seq) < 400

    def test_longest_shorter_than_shortest_fails(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main(
                SIMPLE_FILE,
                "--longest-seq", "10", "--shortest-seq", "100",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )


class TestCLILineLength:
    """--number-lines option."""

    def test_custom_line_length(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--number-lines", "80", "--silent",
            monkeypatch=monkeypatch,
        )
        with open(outfile) as f:
            for line in f:
                if not line.startswith(">"):
                    assert len(line.rstrip("\n")) <= 80

    def test_invalid_number_lines_exits(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main(
                SIMPLE_FILE,
                "--number-lines", "abc", "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )


class TestCLIRandomSubsample:
    """--random-subsample option."""

    def test_subsample(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--random-subsample", "2", "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile)
        assert len(seqs) == 2

    def test_subsample_larger_than_dataset(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        total = len(protfasta.read_fasta(SIMPLE_FILE))
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--random-subsample", str(total + 100), "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile)
        assert len(seqs) == total


class TestCLIPrintStatistics:
    """--print-statistics option."""

    def test_print_statistics(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile, "--print-statistics",
            monkeypatch=monkeypatch,
        )
        assert "[STATS]" in capsys.readouterr().out


class TestCLIRemoveCommaFromHeader:
    """--remove-comma-from-header option."""

    def test_commas_replaced(self, tmp_path, monkeypatch):
        input_file = tmp_path / "comma.fasta"
        input_file.write_text(">protein1, isoform A\nACDEFG\n>protein2, isoform B\nHIKLMN\n")
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            str(input_file), "-o", outfile,
            "--remove-comma-from-header", "--silent",
            monkeypatch=monkeypatch,
        )
        with open(outfile) as f:
            content = f.read()
        assert "," not in content
        assert ";" in content


class TestCLINonUniqueHeader:
    """--non-unique-header option."""

    def test_non_unique_header_allows_duplicates(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            DUPLICATE_RECORD_FILE, "-o", outfile,
            "--non-unique-header", "--duplicate-record", "ignore", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)


class TestCLIInvalidOption:
    """Invalid option values are rejected."""

    def test_bad_duplicate_record_option(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main(
                SIMPLE_FILE,
                "--duplicate-record", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )

    def test_bad_duplicate_sequence_option(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main(
                SIMPLE_FILE,
                "--duplicate-sequence", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )

    def test_bad_invalid_sequence_option(self, monkeypatch):
        with pytest.raises(SystemExit):
            _run_main(
                SIMPLE_FILE,
                "--invalid-sequence", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
