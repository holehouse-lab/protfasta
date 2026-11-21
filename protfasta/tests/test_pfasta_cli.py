"""
Tests for the pfasta command-line interface (protfasta.scripts.pfasta).

Covers:
- Helper functions (stdout, exit_error, validate, validate_int,
  print_statistical_summary)
- CLI argument parsing and end-to-end invocation of main()
"""

from __future__ import annotations

import gzip
import os
import subprocess
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
ALIGNED_VALID_FILE = os.path.join(TEST_DATA_DIR, "aligned_seq_all_valid.fasta")
ALIGNED_CONVERTABLE_FILE = os.path.join(TEST_DATA_DIR, "aligned_seq_all_valid_convertable.fasta")


def _run_main(*cli_args: str, monkeypatch) -> None:
    """Invoke main() after patching sys.argv with the supplied CLI args."""
    monkeypatch.setattr(sys, "argv", ["pfasta", *cli_args])
    main()


def _read_unchecked(filename: str, alignment: bool = False) -> list:
    """Read a FASTA file exactly as it is, with every content check off.

    Used to look at what pfasta wrote (or was given) without the reader
    itself removing, converting or objecting to anything.
    """
    records = protfasta.read_fasta(
        filename,
        expect_unique_header=False,
        duplicate_record_action="ignore",
        duplicate_sequence_action="ignore",
        invalid_sequence_action="ignore",
        alignment=alignment,
        return_list=True,
    )
    assert isinstance(records, list)
    return records


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
        with pytest.raises(SystemExit) as exit_info:
            exit_error("boom")
        assert exit_info.value.code == 1

    def test_exit_error_prints_message(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            exit_error("something broke")
        assert exit_info.value.code == 1
        assert "something broke" in capsys.readouterr().out


class TestValidate:
    """Tests for the validate() helper."""

    def test_valid_option_returned_lowercase(self):
        assert validate("FAIL", ["ignore", "fail", "remove"]) == "fail"

    def test_invalid_option_exits(self):
        with pytest.raises(SystemExit) as exit_info:
            validate("nope", ["ignore", "fail", "remove"])
        assert exit_info.value.code == 1

    def test_mixed_case(self):
        assert validate("Remove", ["ignore", "fail", "remove"]) == "remove"


class TestValidateInt:
    """Tests for the validate_int() helper."""

    def test_valid_int(self):
        assert validate_int("10", 0, "--test") == 10

    def test_value_at_minimum(self):
        assert validate_int("5", 5, "--test") == 5

    def test_below_minimum_exits(self):
        with pytest.raises(SystemExit) as exit_info:
            validate_int("3", 5, "--test")
        assert exit_info.value.code == 1

    def test_non_numeric_exits(self):
        with pytest.raises(SystemExit) as exit_info:
            validate_int("abc", 0, "--test")
        assert exit_info.value.code == 1


class TestPrintStatisticalSummary:
    """Tests for the print_statistical_summary() helper."""

    def test_prints_stats(self, capsys):
        data = [["h1", "ACDE"], ["h2", "ACDEFGHIKLMNPQ"], ["h3", "ACD"]]
        print_statistical_summary(data)
        out = capsys.readouterr().out
        assert "[STATS]" in out
        assert "Total number of sequences" in out

    def test_empty_dataset_does_not_raise(self, capsys):
        print_statistical_summary([])
        assert "[STATS]" in capsys.readouterr().out


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

    def test_no_args_message_has_a_single_error_prefix(self, monkeypatch, capsys):
        # argparse adds 'prog: error:' itself; the message used to carry a
        # second copy
        with pytest.raises(SystemExit) as exit_info:
            _run_main(monkeypatch=monkeypatch)
        assert exit_info.value.code == 2
        err = capsys.readouterr().err
        assert "the following arguments are required: filename" in err
        assert err.count("error:") == 1


class TestCLIMissingFile:
    """Pointing at a non-existent file should error."""

    def test_missing_file(self, monkeypatch):
        with pytest.raises(SystemExit) as exit_info:
            _run_main("no_such_file.fasta", "--silent", monkeypatch=monkeypatch)
        assert exit_info.value.code == 1


class TestCLIUnwritableOutput:
    """An output path that cannot be opened is a clean error, not a traceback."""

    def test_unwritable_output_exits(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "no_such_dir" / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(SIMPLE_FILE, "-o", outfile, "--silent", monkeypatch=monkeypatch)
        assert exc.value.code != 0
        assert "[FATAL ERROR]" in capsys.readouterr().out


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
        # ...and nothing else was written in its place
        assert os.listdir(tmp_path) == []

    def test_output_file_written_without_the_flag(self, tmp_path, monkeypatch):
        # the control for the test above: the same run minus --no-outputfile
        # does write ./output.fasta, so its absence there means something
        monkeypatch.chdir(tmp_path)
        _run_main(SIMPLE_FILE, "--silent", monkeypatch=monkeypatch)
        assert os.listdir(tmp_path) == ["output.fasta"]
        assert _read_unchecked("output.fasta") == _read_unchecked(SIMPLE_FILE)


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
        # three records in, one of them a repeat of another
        written = _read_unchecked(outfile)
        assert len(_read_unchecked(DUPLICATE_RECORD_FILE)) == 3
        assert len(written) == 2
        assert len({(header, seq) for header, seq in written}) == 2

    def test_duplicate_record_fail(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(
                DUPLICATE_RECORD_FILE, "-o", outfile,
                "--duplicate-record", "fail", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code != 0
        assert "[FATAL ERROR]" in capsys.readouterr().out
        assert not os.path.exists(outfile)

    def test_duplicate_record_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            DUPLICATE_RECORD_FILE, "-o", outfile,
            "--duplicate-record", "ignore", "--non-unique-header", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)
        # every record is kept, the repeat included
        assert _read_unchecked(outfile) == _read_unchecked(DUPLICATE_RECORD_FILE)


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
        # three records in, two of them with the same sequence
        written = _read_unchecked(outfile)
        assert len(_read_unchecked(DUPLICATE_SEQ_FILE)) == 3
        assert len(written) == 2
        assert len({seq for _, seq in written}) == 2

    def test_duplicate_sequence_fail(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(
                DUPLICATE_SEQ_FILE, "-o", outfile,
                "--duplicate-sequence", "fail", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code != 0
        assert "[FATAL ERROR]" in capsys.readouterr().out


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
        # written exactly as read, non-standard residues and all
        written = _read_unchecked(outfile)
        assert written == _read_unchecked(NONSTANDARD_FILE)
        assert any("X" in seq for _, seq in written)

    def test_invalid_sequence_fail(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(
                NONSTANDARD_FILE, "-o", outfile,
                "--invalid-sequence", "fail", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code != 0
        out = capsys.readouterr().out
        assert "[FATAL ERROR]" in out
        assert "invalid amino acid" in out

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
        # the same result as the API's convert-ignore, and the convertible
        # residues really are gone
        expected = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action="convert-ignore", return_list=True)
        written = _read_unchecked(outfile)
        assert written == expected
        assert written != _read_unchecked(NONSTANDARD_FILE)
        assert not any(residue in seq for _, seq in written for residue in "BUXZ*- ")

    def test_invalid_sequence_convert_res_ignore(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            NONSTANDARD_FILE, "-o", outfile,
            "--invalid-sequence", "convert-res-ignore", "--silent",
            monkeypatch=monkeypatch,
        )
        assert os.path.exists(outfile)
        # the same result as the API's convert-ignore on an alignment
        expected = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action="convert-ignore",
                                        alignment=True, return_list=True)
        written = _read_unchecked(outfile)
        assert written == expected
        assert not any(residue in seq for _, seq in written for residue in "BUXZ* ")

    # The -res variants exist to handle alignments: '-' must survive as a
    # gap character AND be treated as valid. Previously the CLI stripped
    # '-' from its conversion table but still read the file as a plain
    # (non-alignment) FASTA, so convert-res raised on every gapped
    # sequence and convert-res-remove discarded all of them.

    def test_convert_res_keeps_gaps_on_alignment(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            ALIGNED_CONVERTABLE_FILE, "-o", outfile,
            "--invalid-sequence", "convert-res", "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile, alignment=True)
        assert seqs["Seq1"] == "A-----CDEFGHIKLMNPQRSTVWY"
        assert seqs["Seq2"] == "ACDEFGHIKL-----MNPQRSTVWYN"
        assert seqs["Seq3"] == "ACDEFGHIKLMNPQRSTVWY-----Q"

    def test_convert_res_remove_keeps_aligned_sequences(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            ALIGNED_VALID_FILE, "-o", outfile,
            "--invalid-sequence", "convert-res-remove", "--silent",
            monkeypatch=monkeypatch,
        )
        assert len(protfasta.read_fasta(outfile, alignment=True)) == 3

    def test_convert_all_strips_gaps(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            ALIGNED_VALID_FILE, "-o", outfile,
            "--invalid-sequence", "convert-all", "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile)
        assert seqs["Seq1"] == "ACDEFGHIKLMNPQRSTVWY"
        assert "-" not in "".join(seqs.values())

    def test_convert_all_strips_spaces(self, tmp_path, monkeypatch):
        # Whitespace inside a sequence is convertible via the built-in
        # table (' ' -> ''); the old hand-rolled CLI table omitted it.
        input_file = tmp_path / "spaces.fasta"
        input_file.write_text(">h1\nACD EFG\n>h2\nHIK LMN\n")
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            str(input_file), "-o", outfile,
            "--invalid-sequence", "convert-all", "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile)
        assert seqs == {"h1": "ACDEFG", "h2": "HIKLMN"}


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
        with pytest.raises(SystemExit) as exit_info:
            _run_main(
                SIMPLE_FILE,
                "--longest-seq", "10", "--shortest-seq", "100",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exit_info.value.code == 1

    def test_equal_bounds_fail(self, monkeypatch, capsys):
        # both bounds are exclusive, so equal values can never keep anything
        with pytest.raises(SystemExit) as exc:
            _run_main(
                SIMPLE_FILE,
                "--longest-seq", "100", "--shortest-seq", "100",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code != 0
        assert "must be larger" in capsys.readouterr().out

    def test_shortest_seq_zero_is_honoured(self, tmp_path, monkeypatch):
        # 0 is falsy, but an explicitly-passed 0 must still apply the filter
        # (keeping every sequence of length > 0, i.e. all of them).
        outfile = str(tmp_path / "out.fasta")
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--shortest-seq", "0", "--silent",
            monkeypatch=monkeypatch,
        )
        assert len(protfasta.read_fasta(outfile)) == len(protfasta.read_fasta(SIMPLE_FILE))

    def test_longest_seq_zero_is_honoured(self, tmp_path, monkeypatch):
        # --longest-seq 0 discards everything, so pfasta exits cleanly
        # without writing an output file.
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(
                SIMPLE_FILE, "-o", outfile,
                "--longest-seq", "0", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code == 0
        assert not os.path.exists(outfile)


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
        with pytest.raises(SystemExit) as exit_info:
            _run_main(
                SIMPLE_FILE,
                "--number-lines", "abc", "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exit_info.value.code == 1

    def test_number_lines_zero_is_rejected(self, monkeypatch, capsys):
        # 0 is falsy and used to fall through to the default of 60
        with pytest.raises(SystemExit) as exc:
            _run_main(
                SIMPLE_FILE,
                "--number-lines", "0", "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code != 0
        assert "--number-lines" in capsys.readouterr().out


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

    def test_subsample_is_a_subset_of_input(self, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        full = protfasta.read_fasta(SIMPLE_FILE)
        _run_main(
            SIMPLE_FILE, "-o", outfile,
            "--random-subsample", "3", "--silent",
            monkeypatch=monkeypatch,
        )
        seqs = protfasta.read_fasta(outfile)
        assert len(seqs) == 3
        for header, seq in seqs.items():
            assert full[header] == seq

    def test_subsample_zero_writes_no_file(self, tmp_path, monkeypatch):
        # nothing is left after a subsample of 0, so pfasta says so and
        # exits cleanly rather than writing an empty FASTA file
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(
                SIMPLE_FILE, "-o", outfile,
                "--random-subsample", "0", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exc.value.code == 0
        assert not os.path.exists(outfile)

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
        written = _read_unchecked(outfile)
        assert len(written) == 3
        # a header really is repeated in what was written
        assert len({header for header, _ in written}) == 2

    def test_duplicate_headers_fail_without_the_flag(self, tmp_path, monkeypatch, capsys):
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exit_info:
            _run_main(DUPLICATE_RECORD_FILE, "-o", str(outfile), "--silent", monkeypatch=monkeypatch)
        assert exit_info.value.code == 1
        assert "duplicate header" in capsys.readouterr().out
        assert not outfile.exists()


class TestCLIInvalidOption:
    """Invalid option values are rejected."""

    def test_bad_duplicate_record_option(self, monkeypatch):
        with pytest.raises(SystemExit) as exit_info:
            _run_main(
                SIMPLE_FILE,
                "--duplicate-record", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exit_info.value.code == 1

    def test_bad_duplicate_sequence_option(self, monkeypatch):
        with pytest.raises(SystemExit) as exit_info:
            _run_main(
                SIMPLE_FILE,
                "--duplicate-sequence", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exit_info.value.code == 1

    def test_bad_invalid_sequence_option(self, monkeypatch):
        with pytest.raises(SystemExit) as exit_info:
            _run_main(
                SIMPLE_FILE,
                "--invalid-sequence", "badvalue",
                "--no-outputfile", "--silent",
                monkeypatch=monkeypatch,
            )
        assert exit_info.value.code == 1


class TestCLIEmptySequence:
    """--empty-sequence controls headers with no sequence after them."""

    @pytest.fixture
    def empty_file(self, tmp_path):
        f = tmp_path / "empty_records.fasta"
        f.write_text(">good\nACD\n>empty_one\n>also_good\nEFG\n")
        return str(f)

    def test_fails_by_default(self, empty_file, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exc:
            _run_main(empty_file, "-o", outfile, "--silent", monkeypatch=monkeypatch)
        assert exc.value.code != 0
        out = capsys.readouterr().out
        assert "[FATAL ERROR]" in out
        assert "empty_one" in out

    def test_remove(self, empty_file, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(empty_file, "-o", outfile, "--empty-sequence", "remove", "--silent", monkeypatch=monkeypatch)
        assert protfasta.read_fasta(outfile) == {"good": "ACD", "also_good": "EFG"}

    def test_ignore_not_offered(self, empty_file, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        with pytest.raises(SystemExit) as exit_info:
            _run_main(empty_file, "-o", outfile, "--empty-sequence", "ignore", "--silent", monkeypatch=monkeypatch)
        assert exit_info.value.code == 1


class TestCLILeadingContent:
    """--leading-content controls content before the first header."""

    @pytest.fixture
    def leading_file(self, tmp_path):
        f = tmp_path / "leading_content.fasta"
        f.write_text("junk\n>a\nACD\n>b\nEFG\n")
        return str(f)

    def test_fails_by_default(self, leading_file, tmp_path, monkeypatch, capsys):
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exc:
            _run_main(leading_file, "-o", str(outfile), "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        out = capsys.readouterr().out
        assert "[FATAL ERROR]" in out
        assert "before the first header" in out
        assert not outfile.exists()

    def test_compressed_input_is_an_error_not_zero_sequences(self, tmp_path, monkeypatch, capsys):
        # this used to report "0 sequences remain after filtering" and exit 0
        infile = tmp_path / "six.fasta.gz"
        infile.write_bytes(gzip.compress(b">a\nACD\n>b\nEFG\n"))
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exc:
            _run_main(str(infile), "-o", str(outfile), monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "before the first header" in capsys.readouterr().out
        assert not outfile.exists()

    def test_remove(self, leading_file, tmp_path, monkeypatch):
        outfile = str(tmp_path / "out.fasta")
        _run_main(leading_file, "-o", outfile, "--leading-content", "remove", "--silent", monkeypatch=monkeypatch)
        assert protfasta.read_fasta(outfile) == {"a": "ACD", "b": "EFG"}

    def test_invalid_value_rejected(self, leading_file, tmp_path, monkeypatch):
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exc:
            _run_main(leading_file, "-o", str(outfile), "--leading-content", "ignore", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert not outfile.exists()


def _write_sequences_of_lengths(path, lengths) -> str:
    """Write a FASTA file with one poly-alanine record per requested length."""
    with open(path, "w") as fh:
        for length in lengths:
            fh.write(">len%i\n%s\n" % (length, "A" * length))
    return str(path)


class TestCLIRefusesToOverwriteInput:
    """pfasta must never write its output over its own input."""

    def test_explicit_output_equal_to_input(self, tmp_path, monkeypatch, capsys):
        infile = _write_sequences_of_lengths(tmp_path / "same.fasta", [1, 2, 3, 4, 5, 6])
        before = (tmp_path / "same.fasta").read_bytes()
        with pytest.raises(SystemExit) as exc:
            _run_main(infile, "-o", infile, "--shortest-seq", "4", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "is the input file" in capsys.readouterr().out
        assert (tmp_path / "same.fasta").read_bytes() == before

    def test_default_output_name_equal_to_input(self, tmp_path, monkeypatch):
        # with no -o the output is ./output.fasta, which here is the input
        _write_sequences_of_lengths(tmp_path / "output.fasta", [1, 2, 3, 4, 5, 6])
        before = (tmp_path / "output.fasta").read_bytes()
        monkeypatch.chdir(tmp_path)
        with pytest.raises(SystemExit) as exc:
            _run_main("output.fasta", "--longest-seq", "3", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert (tmp_path / "output.fasta").read_bytes() == before

    def test_symlink_to_input(self, tmp_path, monkeypatch):
        infile = _write_sequences_of_lengths(tmp_path / "real.fasta", [1, 2, 3])
        link = tmp_path / "link.fasta"
        link.symlink_to(tmp_path / "real.fasta")
        before = (tmp_path / "real.fasta").read_bytes()
        with pytest.raises(SystemExit) as exc:
            _run_main(infile, "-o", str(link), "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert (tmp_path / "real.fasta").read_bytes() == before

    def test_no_outputfile_run_on_output_fasta_is_allowed(self, tmp_path, monkeypatch, capsys):
        # nothing is written, so there is nothing to protect
        _write_sequences_of_lengths(tmp_path / "output.fasta", [1, 2, 3])
        monkeypatch.chdir(tmp_path)
        _run_main("output.fasta", "--no-outputfile", "--print-statistics", monkeypatch=monkeypatch)
        assert "Total number of sequences : 3" in capsys.readouterr().out


class TestCLIStatisticsValues:
    """--print-statistics reports interpolated, unbiased quantiles."""

    def _stats(self, tmp_path, monkeypatch, capsys, lengths) -> str:
        infile = _write_sequences_of_lengths(tmp_path / "lengths.fasta", lengths)
        _run_main(infile, "--no-outputfile", "--print-statistics", monkeypatch=monkeypatch)
        return capsys.readouterr().out

    def test_median_of_two(self, tmp_path, monkeypatch, capsys):
        out = self._stats(tmp_path, monkeypatch, capsys, [10, 20])
        assert "Median length             : 15.0" in out
        # two sequences are enough to interpolate the quartiles as well
        assert "25th percentile length    : 12.5" in out
        assert "75th percentile length    : 17.5" in out

    def test_quartiles_of_four(self, tmp_path, monkeypatch, capsys):
        out = self._stats(tmp_path, monkeypatch, capsys, [10, 20, 30, 40])
        assert "25th percentile length    : 17.5" in out
        assert "Median length             : 25.0" in out
        assert "75th percentile length    : 32.5" in out
        assert "Longest sequence          : 40" in out
        assert "Shortest sequence         : 10" in out

    def test_single_sequence(self, tmp_path, monkeypatch, capsys):
        out = self._stats(tmp_path, monkeypatch, capsys, [7])
        assert "25th percentile length    : 7.0" in out
        assert "Median length             : 7.0" in out
        assert "75th percentile length    : 7.0" in out


class TestCLILengthFilterBoundaries:
    """Both length bounds are exclusive, tested at the boundary itself."""

    def test_bounds_are_exclusive(self, tmp_path, monkeypatch):
        infile = _write_sequences_of_lengths(tmp_path / "lengths.fasta", [1, 2, 3, 4, 5, 6])
        outfile = str(tmp_path / "out.fasta")
        _run_main(infile, "-o", outfile, "--shortest-seq", "2", "--longest-seq", "5", "--silent", monkeypatch=monkeypatch)
        assert sorted(protfasta.read_fasta(outfile)) == ["len3", "len4"]

    def test_adjacent_bounds_fail(self, tmp_path, monkeypatch, capsys):
        # nothing is longer than 3 and shorter than 4
        infile = _write_sequences_of_lengths(tmp_path / "lengths.fasta", [3, 4])
        with pytest.raises(SystemExit) as exc:
            _run_main(infile, "--shortest-seq", "3", "--longest-seq", "4", "--no-outputfile", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "must be larger" in capsys.readouterr().out

    def test_bounds_two_apart_are_accepted(self, tmp_path, monkeypatch):
        infile = _write_sequences_of_lengths(tmp_path / "lengths.fasta", [3, 4, 5])
        outfile = str(tmp_path / "out.fasta")
        _run_main(infile, "-o", outfile, "--shortest-seq", "3", "--longest-seq", "5", "--silent", monkeypatch=monkeypatch)
        assert list(protfasta.read_fasta(outfile)) == ["len4"]


class TestCLIRandomSeed:
    """--random-seed makes --random-subsample reproducible."""

    def _subsample(self, tmp_path, monkeypatch, name, *extra) -> list:
        outfile = str(tmp_path / name)
        _run_main(SIMPLE_FILE, "-o", outfile, "--random-subsample", "4", "--silent", *extra, monkeypatch=monkeypatch)
        records = protfasta.read_fasta(outfile, return_list=True)
        assert isinstance(records, list)
        return records

    def test_same_seed_same_subsample(self, tmp_path, monkeypatch):
        first = self._subsample(tmp_path, monkeypatch, "first.fasta", "--random-seed", "42")
        second = self._subsample(tmp_path, monkeypatch, "second.fasta", "--random-seed", "42")
        assert len(first) == 4
        assert first == second

    def test_seed_selects_the_subsample(self, tmp_path, monkeypatch):
        # nine sequences give 3024 ordered sub-samples of four, so twenty
        # different seeds all drawing the same one would mean the seed is
        # being ignored
        draws = [self._subsample(tmp_path, monkeypatch, "seed%i.fasta" % seed, "--random-seed", str(seed)) for seed in range(20)]
        assert any(draw != draws[0] for draw in draws)

    @pytest.mark.parametrize("bad_seed", ["abc", "-1", "1.5"])
    def test_invalid_seed_rejected(self, tmp_path, monkeypatch, bad_seed):
        with pytest.raises(SystemExit) as exc:
            self._subsample(tmp_path, monkeypatch, "out.fasta", "--random-seed", bad_seed)
        assert exc.value.code == 1


class TestCLIOptionValues:
    """Options do what they say, and an empty value is not 'use the default'."""

    @pytest.mark.parametrize("flag", ["--duplicate-record", "--duplicate-sequence", "--empty-sequence", "--leading-content", "--invalid-sequence"])
    def test_empty_string_value_rejected(self, tmp_path, monkeypatch, flag):
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exc:
            _run_main(SIMPLE_FILE, "-o", str(outfile), flag, "", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert not outfile.exists()

    def test_empty_output_name_rejected(self, monkeypatch, capsys):
        with pytest.raises(SystemExit) as exc:
            _run_main(SIMPLE_FILE, "-o", "", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "[FATAL ERROR]" in capsys.readouterr().out

    def test_invalid_residues_fail_by_default(self, tmp_path, monkeypatch, capsys):
        infile = tmp_path / "nonstandard.fasta"
        infile.write_text(">a\nACX\n")
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exc:
            _run_main(str(infile), "-o", str(outfile), "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "invalid amino acid: X" in capsys.readouterr().out
        assert not outfile.exists()

    def test_duplicate_sequence_remove_removes(self, tmp_path, monkeypatch):
        infile = tmp_path / "duplicates.fasta"
        infile.write_text(">a\nACD\n>b\nACD\n>c\nEFG\n")
        outfile = str(tmp_path / "out.fasta")
        _run_main(str(infile), "-o", outfile, "--duplicate-sequence", "remove", "--silent", monkeypatch=monkeypatch)
        assert protfasta.read_fasta(outfile) == {"a": "ACD", "c": "EFG"}

    def test_sequence_emptied_by_conversion_fails(self, tmp_path, monkeypatch, capsys):
        infile = tmp_path / "stars.fasta"
        infile.write_text(">a\n***\n>b\nACD\n")
        with pytest.raises(SystemExit) as exc:
            _run_main(str(infile), "--invalid-sequence", "convert-all", "--no-outputfile", "--silent", monkeypatch=monkeypatch)
        assert exc.value.code == 1
        assert "left no sequence" in capsys.readouterr().out


class TestCLIAsModule:
    """`python -m protfasta.scripts.pfasta` runs the tool (it once did nothing)."""

    def test_version_via_python_dash_m(self):
        completed = subprocess.run(
            [sys.executable, "-m", "protfasta.scripts.pfasta", "--version"],
            capture_output=True, text=True,
        )
        assert completed.returncode == 0
        assert completed.stdout.strip() == VERSION_MAJ


class TestCLIDefaults:
    """What pfasta does when an option is not given."""

    def test_gapped_alignment_fails_by_default(self, tmp_path, monkeypatch, capsys):
        # a '-' is only acceptable when a convert-res mode asks for alignment handling
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exit_info:
            _run_main(ALIGNED_VALID_FILE, "-o", str(outfile), "--silent", monkeypatch=monkeypatch)
        assert exit_info.value.code == 1
        assert "invalid amino acid: -" in capsys.readouterr().out
        assert not outfile.exists()

    def test_default_line_length_is_60(self, tmp_path, monkeypatch):
        outfile = tmp_path / "out.fasta"
        _run_main(SIMPLE_FILE, "-o", str(outfile), "--silent", monkeypatch=monkeypatch)
        sequence_lines = [line for line in outfile.read_text().split("\n") if line and not line.startswith(">")]
        assert max(len(line) for line in sequence_lines) == 60

    @pytest.mark.parametrize("number_lines", ["5", "6", "100"])
    def test_number_lines_sets_the_line_length(self, tmp_path, monkeypatch, number_lines):
        outfile = tmp_path / "out.fasta"
        _run_main(SIMPLE_FILE, "-o", str(outfile), "--number-lines", number_lines, "--silent", monkeypatch=monkeypatch)
        sequence_lines = [line for line in outfile.read_text().split("\n") if line and not line.startswith(">")]
        assert max(len(line) for line in sequence_lines) == int(number_lines)
        assert _read_unchecked(str(outfile)) == _read_unchecked(SIMPLE_FILE)

    def test_number_lines_below_five_rejected(self, tmp_path, monkeypatch, capsys):
        outfile = tmp_path / "out.fasta"
        with pytest.raises(SystemExit) as exit_info:
            _run_main(SIMPLE_FILE, "-o", str(outfile), "--number-lines", "4", "--silent", monkeypatch=monkeypatch)
        assert exit_info.value.code == 1
        assert "--number-lines must be a numerical value >= 5" in capsys.readouterr().out
        assert not outfile.exists()

    def test_no_statistics_unless_asked(self, tmp_path, monkeypatch, capsys):
        _run_main(SIMPLE_FILE, "-o", str(tmp_path / "out.fasta"), monkeypatch=monkeypatch)
        assert "[STATS]" not in capsys.readouterr().out

    def test_silent_wins_over_print_statistics(self, tmp_path, monkeypatch, capsys):
        _run_main(SIMPLE_FILE, "-o", str(tmp_path / "out.fasta"), "--print-statistics", "--silent", monkeypatch=monkeypatch)
        assert capsys.readouterr().out == ""

    def test_progress_is_reported_unless_silent(self, tmp_path, monkeypatch, capsys):
        outfile = str(tmp_path / "out.fasta")
        _run_main(SIMPLE_FILE, "-o", outfile, monkeypatch=monkeypatch)
        out = capsys.readouterr().out
        assert "Reading in the file %s" % SIMPLE_FILE in out
        assert "Parsed file to recover 9 sequences" in out
        assert "Writing new sequence file [%s]" % outfile in out

    def test_statistics_report_the_real_extremes(self, tmp_path, monkeypatch, capsys):
        infile = _write_sequences_of_lengths(tmp_path / "lengths.fasta", [3, 9, 4, 7, 5])
        _run_main(infile, "--no-outputfile", "--print-statistics", monkeypatch=monkeypatch)
        out = capsys.readouterr().out
        assert "Total number of sequences : 5" in out
        assert "25th percentile length    : 4.0" in out
        assert "Median length             : 5.0" in out
        assert "75th percentile length    : 7.0" in out
        assert "Longest sequence          : 9" in out
        assert "Shortest sequence         : 3" in out


class TestCLISubsampleMessages:
    """--random-subsample says what it is doing, including when it can't."""

    def _run(self, tmp_path, monkeypatch, capsys, n_requested: int) -> tuple:
        outfile = str(tmp_path / "out.fasta")
        _run_main(SIMPLE_FILE, "-o", outfile, "--random-subsample", str(n_requested), "--random-seed", "1", monkeypatch=monkeypatch)
        return capsys.readouterr().out, _read_unchecked(outfile)

    def test_fewer_than_available(self, tmp_path, monkeypatch, capsys):
        out, written = self._run(tmp_path, monkeypatch, capsys, 4)
        assert "Subsampling 4 sequences from the complete dataset (9)" in out
        assert len(written) == 4
        # a sub-sample is drawn without replacement from the input
        original = _read_unchecked(SIMPLE_FILE)
        assert all(record in original for record in written)
        assert len({header for header, _ in written}) == 4

    def test_exactly_as_many_as_available(self, tmp_path, monkeypatch, capsys):
        out, written = self._run(tmp_path, monkeypatch, capsys, 9)
        assert "Subsampling 9 sequences from the complete dataset (9)" in out
        assert "Cannot subsample" not in out
        assert sorted(written) == sorted(_read_unchecked(SIMPLE_FILE))

    def test_more_than_available(self, tmp_path, monkeypatch, capsys):
        out, written = self._run(tmp_path, monkeypatch, capsys, 10)
        assert "Cannot subsample" in out
        assert "Subsampling" not in out
        assert sorted(written) == sorted(_read_unchecked(SIMPLE_FILE))
