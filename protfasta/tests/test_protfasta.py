"""
Systematic test suite for the protfasta package.

Organized into test classes by functional area:
- TestImport: Package import verification
- TestConfigs: Configuration constants
- TestCheckSequenceIsValid: Sequence validation utility
- TestConvertToValid: Sequence conversion utility
- TestBuildCustomDictionary: Dictionary building utility
- TestRemoveInvalidSequences: Invalid sequence removal
- TestFailOnInvalidSequences: Invalid sequence failure
- TestConvertInvalidSequences: Invalid sequence conversion
- TestDuplicateUtilities: Duplicate detection/removal utilities
- TestConvertListToDictionary: List to dict conversion
- TestCheckInputs: Input parameter validation
- TestInternalParseFasta: Low-level FASTA parsing
- TestReadFastaBasic: Basic read_fasta usage
- TestReadFastaExpectUniqueHeader: expect_unique_header parameter
- TestReadFastaHeaderParser: header_parser parameter
- TestReadFastaDuplicateRecords: duplicate_record_action parameter
- TestReadFastaDuplicateSequences: duplicate_sequence_action parameter
- TestReadFastaInvalidSequences: invalid_sequence_action parameter
- TestReadFastaAlignment: alignment parameter
- TestReadFastaReturnList: return_list parameter
- TestReadFastaOutputFile: output_filename parameter
- TestReadFastaCorrectionDictionary: correction_dictionary parameter
- TestWriteFasta: write_fasta function
- TestEndToEndCombinations: Combined parameter interactions
"""

import protfasta
from protfasta.protfasta_exceptions import ProtfastaException
from protfasta import _configs
from protfasta import utilities as _utilities
from protfasta import io as _io
import pytest
import sys
import os

from pathlib import Path

HERE = Path(__file__).parent
TEST_DATA_DIR = protfasta._get_data('test_data')

# Known reference data from testset_1.fasta
WASL_HEADER = 'sp|O00401|WASL_HUMAN Neural Wiskott-Aldrich syndrome protein OS=Homo sapiens OX=9606 GN=WASL PE=1 SV=2'
WASL_SEQ = (
    'MSSVQQQPPPPRRVTNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKK'
    'CSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALN'
    'FANEEEAKKFRKAVTDLLGRRQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVN'
    'NISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGWDPNTGFDLNNLDPELKNLFDMCG'
    'ISEAQLKDRETSKVIYDFIEKTGGVEAVKNELRRQAPPPPPPSRGGPPPPPPPPHNSGPP'
    'PPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPP'
    'PPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLK'
    'KVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSK'
    'AIHSSDEDEDEDDEEDFEDDDEWED'
)

# File paths
SIMPLE_FILE = os.path.join(TEST_DATA_DIR, 'testset_1.fasta')
DUPLICATE_RECORD_FILE = os.path.join(TEST_DATA_DIR, 'testset_duplicate.fasta')
DUPLICATE_SEQ_FILE = os.path.join(TEST_DATA_DIR, 'testset_duplicate_seqs.fasta')
BADCHAR_FILE = os.path.join(TEST_DATA_DIR, 'test_data_with_bad_chars.fa')
NONSTANDARD_FILE = os.path.join(TEST_DATA_DIR, 'test_data_with_nonstandard_chars.fa')
FIXABLE_INVALID_FILE = os.path.join(TEST_DATA_DIR, 'fixable_invalid.fasta')
UNFIXABLE_INVALID_FILE = os.path.join(TEST_DATA_DIR, 'unfixable_invalid.fasta')
ALIGNED_VALID_FILE = os.path.join(TEST_DATA_DIR, 'aligned_seq_all_valid.fasta')
ALIGNED_CONVERTABLE_FILE = os.path.join(TEST_DATA_DIR, 'aligned_seq_all_valid_convertable.fasta')
ALIGNED_UNCONVERTABLE_FILE = os.path.join(TEST_DATA_DIR, 'aligned_seq_all_valid_unconvertable.fasta')


# ---------------------------------------------------------------------------
# TestImport
# ---------------------------------------------------------------------------
class TestImport:
    """Verify that the package can be imported and exposes the expected API."""

    def test_module_imported(self):
        assert "protfasta" in sys.modules

    def test_read_fasta_callable(self):
        assert callable(protfasta.read_fasta)

    def test_write_fasta_callable(self):
        assert callable(protfasta.write_fasta)

    def test_version_exists(self):
        assert hasattr(protfasta, '__version__')
        assert isinstance(protfasta.__version__, str)

    def test_get_data_returns_valid_path(self):
        path = protfasta._get_data('test_data')
        assert os.path.isdir(path)

    def test_exception_class_importable(self):
        assert issubclass(ProtfastaException, Exception)


# ---------------------------------------------------------------------------
# TestConfigs
# ---------------------------------------------------------------------------
class TestConfigs:
    """Verify configuration constants are correctly defined."""

    def test_standard_aas_count(self):
        assert len(_configs.STANDARD_AAS) == 20

    def test_standard_aas_are_uppercase_single_chars(self):
        for aa in _configs.STANDARD_AAS:
            assert len(aa) == 1
            assert aa.isupper()

    def test_standard_aas_with_gap_includes_dash(self):
        assert '-' in _configs.STANDARD_AAS_WITH_GAP
        assert len(_configs.STANDARD_AAS_WITH_GAP) == 21

    def test_standard_conversion_keys(self):
        expected_keys = {'B', 'U', 'X', 'Z', '*', '-', ' '}
        assert set(_configs.STANDARD_CONVERSION.keys()) == expected_keys

    def test_standard_conversion_values(self):
        assert _configs.STANDARD_CONVERSION['B'] == 'N'
        assert _configs.STANDARD_CONVERSION['U'] == 'C'
        assert _configs.STANDARD_CONVERSION['X'] == 'G'
        assert _configs.STANDARD_CONVERSION['Z'] == 'Q'
        assert _configs.STANDARD_CONVERSION['*'] == ''
        assert _configs.STANDARD_CONVERSION['-'] == ''
        assert _configs.STANDARD_CONVERSION[' '] == ''

    def test_standard_conversion_with_gap_preserves_dash(self):
        # STANDARD_CONVERSION_WITH_GAP should NOT have '-' as a key
        assert '-' not in _configs.STANDARD_CONVERSION_WITH_GAP

    def test_standard_conversion_with_gap_keys(self):
        expected_keys = {'B', 'U', 'X', 'Z', ' ', '*'}
        assert set(_configs.STANDARD_CONVERSION_WITH_GAP.keys()) == expected_keys


# ---------------------------------------------------------------------------
# TestCheckSequenceIsValid
# ---------------------------------------------------------------------------
class TestCheckSequenceIsValid:
    """Tests for utilities.check_sequence_is_valid."""

    def test_all_standard_aas_valid(self):
        seq = ''.join(_configs.STANDARD_AAS)
        valid, info = _utilities.check_sequence_is_valid(seq)
        assert valid is True
        assert info == 0

    def test_single_valid_residue(self):
        for aa in _configs.STANDARD_AAS:
            valid, info = _utilities.check_sequence_is_valid(aa)
            assert valid is True

    def test_invalid_residue_detected(self):
        valid, info = _utilities.check_sequence_is_valid('ACDEFX')
        assert valid is False
        assert info == 'X'

    def test_nonstandard_B(self):
        valid, info = _utilities.check_sequence_is_valid('ACDB')
        assert valid is False
        assert info == 'B'

    def test_nonstandard_U(self):
        valid, info = _utilities.check_sequence_is_valid('ACDU')
        assert valid is False
        assert info == 'U'

    def test_asterisk_invalid(self):
        valid, info = _utilities.check_sequence_is_valid('ACD*')
        assert valid is False
        assert info == '*'

    def test_dash_invalid_without_alignment(self):
        valid, info = _utilities.check_sequence_is_valid('A--CD')
        assert valid is False
        assert info == '-'

    def test_dash_valid_with_alignment(self):
        valid, info = _utilities.check_sequence_is_valid('A--CD', alignment=True)
        assert valid is True
        assert info == 0

    def test_empty_sequence_is_valid(self):
        valid, info = _utilities.check_sequence_is_valid('')
        assert valid is True

    def test_completely_invalid_char(self):
        valid, info = _utilities.check_sequence_is_valid('ACD@')
        assert valid is False
        assert info == '@'

    def test_lowercase_invalid(self):
        # sequences should be uppercased before validation by the parser
        valid, info = _utilities.check_sequence_is_valid('acde')
        assert valid is False


# ---------------------------------------------------------------------------
# TestConvertToValid
# ---------------------------------------------------------------------------
class TestConvertToValid:
    """Tests for utilities.convert_to_valid."""

    def test_no_change_for_valid_sequence(self):
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        assert _utilities.convert_to_valid(seq) == seq

    def test_B_converted_to_N(self):
        assert _utilities.convert_to_valid('ACDB') == 'ACDN'

    def test_U_converted_to_C(self):
        assert _utilities.convert_to_valid('ACDU') == 'ACDC'

    def test_X_converted_to_G(self):
        assert _utilities.convert_to_valid('ACDX') == 'ACDG'

    def test_Z_converted_to_Q(self):
        assert _utilities.convert_to_valid('ACDZ') == 'ACDQ'

    def test_asterisk_removed(self):
        assert _utilities.convert_to_valid('ACD*') == 'ACD'

    def test_dash_removed_without_alignment(self):
        assert _utilities.convert_to_valid('A--CD') == 'ACD'

    def test_dash_preserved_with_alignment(self):
        assert _utilities.convert_to_valid('A--CD', alignment=True) == 'A--CD'

    def test_space_removed(self):
        assert _utilities.convert_to_valid('A C D') == 'ACD'

    def test_multiple_conversions(self):
        assert _utilities.convert_to_valid('BUXZ') == 'NCGQ'

    def test_custom_dictionary(self):
        cd = {'.': 'A', '-': 'C'}
        assert _utilities.convert_to_valid('G.H-I', correction_dictionary=cd) == 'GAHCI'

    def test_custom_dictionary_overrides_defaults(self):
        # Custom dictionary should completely replace defaults
        cd = {'X': 'A'}
        result = _utilities.convert_to_valid('AXB', correction_dictionary=cd)
        # Only X is converted because custom dict replaces the default
        assert result == 'AAB'


# ---------------------------------------------------------------------------
# TestBuildCustomDictionary
# ---------------------------------------------------------------------------
class TestBuildCustomDictionary:
    """Tests for utilities.build_custom_dictionary."""

    def test_empty_additional_dict(self):
        result = _utilities.build_custom_dictionary({})
        assert result == _configs.STANDARD_CONVERSION

    def test_additional_entries_merged(self):
        result = _utilities.build_custom_dictionary({'.': 'A'})
        assert result['.'] == 'A'
        # standard entries still present
        assert result['B'] == 'N'

    def test_override_standard_entry(self):
        result = _utilities.build_custom_dictionary({'B': 'G'})
        assert result['B'] == 'G'  # overridden, not 'N'

    def test_all_standard_keys_retained(self):
        result = _utilities.build_custom_dictionary({'Q': 'Z'})
        for key in _configs.STANDARD_CONVERSION:
            assert key in result


# ---------------------------------------------------------------------------
# TestRemoveInvalidSequences
# ---------------------------------------------------------------------------
class TestRemoveInvalidSequences:
    """Tests for utilities.remove_invalid_sequences."""

    def test_all_valid_kept(self):
        data = [['h1', 'ACDEF'], ['h2', 'GHIKL']]
        result = _utilities.remove_invalid_sequences(data)
        assert len(result) == 2

    def test_invalid_removed(self):
        data = [['h1', 'ACDEF'], ['h2', 'GHXKL']]
        result = _utilities.remove_invalid_sequences(data)
        assert len(result) == 1
        assert result[0][0] == 'h1'

    def test_all_invalid_gives_empty(self):
        data = [['h1', 'XBUZ'], ['h2', 'X*-']]
        result = _utilities.remove_invalid_sequences(data)
        assert len(result) == 0

    def test_alignment_mode_keeps_dashes(self):
        data = [['h1', 'A--CD'], ['h2', 'GHXKL']]
        result = _utilities.remove_invalid_sequences(data, alignment=True)
        assert len(result) == 1
        assert result[0][0] == 'h1'


# ---------------------------------------------------------------------------
# TestFailOnInvalidSequences
# ---------------------------------------------------------------------------
class TestFailOnInvalidSequences:
    """Tests for utilities.fail_on_invalid_sequences."""

    def test_valid_sequences_no_error(self):
        data = [['h1', 'ACDEF'], ['h2', 'GHIKL']]
        _utilities.fail_on_invalid_sequences(data)  # should not raise

    def test_invalid_sequence_raises(self):
        data = [['h1', 'ACXEF']]
        with pytest.raises(ProtfastaException):
            _utilities.fail_on_invalid_sequences(data)

    def test_dash_raises_without_alignment(self):
        data = [['h1', 'A--CD']]
        with pytest.raises(ProtfastaException):
            _utilities.fail_on_invalid_sequences(data)

    def test_dash_ok_with_alignment(self):
        data = [['h1', 'A--CD']]
        _utilities.fail_on_invalid_sequences(data, alignment=True)


# ---------------------------------------------------------------------------
# TestConvertInvalidSequences
# ---------------------------------------------------------------------------
class TestConvertInvalidSequences:
    """Tests for utilities.convert_invalid_sequences."""

    def test_converts_nonstandard(self):
        data = [['h1', 'ACBX']]
        result, count = _utilities.convert_invalid_sequences(data)
        assert result[0][1] == 'ACNG'
        assert count == 1

    def test_no_conversion_needed(self):
        data = [['h1', 'ACDEF']]
        result, count = _utilities.convert_invalid_sequences(data)
        assert result[0][1] == 'ACDEF'
        assert count == 0

    def test_count_reflects_modified_sequences(self):
        data = [['h1', 'ACBX'], ['h2', 'ACDEF'], ['h3', 'UZW']]
        result, count = _utilities.convert_invalid_sequences(data)
        assert count == 2

    def test_custom_dictionary(self):
        data = [['h1', 'AC.F']]
        result, count = _utilities.convert_invalid_sequences(data, correction_dictionary={'.': 'D'})
        assert result[0][1] == 'ACDF'
        assert count == 1


# ---------------------------------------------------------------------------
# TestDuplicateUtilities
# ---------------------------------------------------------------------------
class TestDuplicateUtilities:
    """Tests for duplicate detection/removal utilities."""

    def test_fail_on_duplicates_no_duplicates(self):
        data = [['h1', 'ACDEF'], ['h2', 'GHIKL']]
        _utilities.fail_on_duplicates(data)  # should not raise

    def test_fail_on_duplicates_with_duplicate(self):
        data = [['h1', 'ACDEF'], ['h1', 'ACDEF']]
        with pytest.raises(ProtfastaException):
            _utilities.fail_on_duplicates(data)

    def test_fail_on_duplicates_same_header_different_seq(self):
        # Same header but different sequence - not a duplicate record
        data = [['h1', 'ACDEF'], ['h1', 'GHIKL']]
        _utilities.fail_on_duplicates(data)  # should not raise

    def test_remove_duplicates(self):
        data = [['h1', 'ACDEF'], ['h1', 'ACDEF'], ['h2', 'GHIKL']]
        result = _utilities.remove_duplicates(data)
        assert len(result) == 2

    def test_remove_duplicates_keeps_first(self):
        data = [['h1', 'SEQ1'], ['h1', 'SEQ1'], ['h1', 'SEQ2']]
        result = _utilities.remove_duplicates(data)
        assert len(result) == 2
        assert result[0][1] == 'SEQ1'
        assert result[1][1] == 'SEQ2'

    def test_fail_on_duplicate_sequences_no_duplicates(self):
        data = [['h1', 'ACDEF'], ['h2', 'GHIKL']]
        _utilities.fail_on_duplicate_sequences(data)  # should not raise

    def test_fail_on_duplicate_sequences_with_duplicate(self):
        data = [['h1', 'ACDEF'], ['h2', 'ACDEF']]
        with pytest.raises(ProtfastaException):
            _utilities.fail_on_duplicate_sequences(data)

    def test_remove_duplicate_sequences(self):
        data = [['h1', 'ACDEF'], ['h2', 'ACDEF'], ['h3', 'GHIKL']]
        result = _utilities.remove_duplicate_sequences(data)
        assert len(result) == 2
        assert result[0][0] == 'h1'
        assert result[1][0] == 'h3'

    def test_remove_duplicate_sequences_keeps_first(self):
        data = [['h1', 'SAME'], ['h2', 'SAME'], ['h3', 'SAME']]
        result = _utilities.remove_duplicate_sequences(data)
        assert len(result) == 1
        assert result[0][0] == 'h1'


# ---------------------------------------------------------------------------
# TestConvertListToDictionary
# ---------------------------------------------------------------------------
class TestConvertListToDictionary:
    """Tests for utilities.convert_list_to_dictionary."""

    def test_basic_conversion(self):
        raw = [['h1', 'ACDEF'], ['h2', 'GHIKL']]
        result = _utilities.convert_list_to_dictionary(raw)
        assert result == {'h1': 'ACDEF', 'h2': 'GHIKL'}

    def test_empty_list(self):
        result = _utilities.convert_list_to_dictionary([])
        assert result == {}

    def test_duplicate_headers_last_wins(self):
        raw = [['h1', 'SEQ1'], ['h1', 'SEQ2']]
        result = _utilities.convert_list_to_dictionary(raw)
        assert result['h1'] == 'SEQ2'

    def test_verbose_mode(self, capsys):
        raw = [['h1', 'ACDEF']]
        _utilities.convert_list_to_dictionary(raw, verbose=True)
        captured = capsys.readouterr()
        assert 'uniquely added' in captured.out

    def test_verbose_with_duplicates(self, capsys):
        raw = [['h1', 'SEQ1'], ['h1', 'SEQ2']]
        _utilities.convert_list_to_dictionary(raw, verbose=True)
        captured = capsys.readouterr()
        assert 'Overwriting' in captured.out


# ---------------------------------------------------------------------------
# TestCheckInputs
# ---------------------------------------------------------------------------
class TestCheckInputs:
    """Tests for io.check_inputs validation."""

    def _call_check_inputs(self, **kwargs):
        """Helper with valid defaults for all parameters."""
        defaults = dict(
            expect_unique_header=True,
            header_parser=None,
            check_header_parser=True,
            duplicate_record_action='fail',
            duplicate_sequence_action='ignore',
            invalid_sequence_action='fail',
            alignment=False,
            return_list=False,
            output_filename=None,
            verbose=False,
            correction_dictionary=None,
        )
        defaults.update(kwargs)
        _io.check_inputs(**defaults)

    def test_valid_defaults_pass(self):
        self._call_check_inputs()

    def test_expect_unique_header_non_bool_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(expect_unique_header='yes')
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(expect_unique_header=1)
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(expect_unique_header=None)

    def test_header_parser_non_callable_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(header_parser='not_a_function')

    def test_header_parser_wrong_signature_fails(self):
        def bad_func():
            return "hello"
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(header_parser=bad_func)

    def test_header_parser_returns_non_string_fails(self):
        def bad_func(s):
            return 42
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(header_parser=bad_func)

    def test_header_parser_valid_passes(self):
        def good_func(s):
            return s.upper()
        self._call_check_inputs(header_parser=good_func)

    def test_header_parser_check_disabled(self):
        # A bad function should be accepted when check_header_parser=False
        def bad_func():
            return "hello"
        self._call_check_inputs(header_parser=bad_func, check_header_parser=False)

    def test_duplicate_record_action_invalid_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(duplicate_record_action='invalid')
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(duplicate_record_action='FAIL')

    def test_duplicate_record_action_valid_values(self):
        for val in ['ignore', 'fail', 'remove']:
            if val == 'ignore':
                self._call_check_inputs(duplicate_record_action=val, expect_unique_header=False)
            else:
                self._call_check_inputs(duplicate_record_action=val)

    def test_duplicate_sequence_action_invalid_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(duplicate_sequence_action='invalid')

    def test_duplicate_sequence_action_valid_values(self):
        for val in ['ignore', 'fail', 'remove']:
            self._call_check_inputs(duplicate_sequence_action=val)

    def test_invalid_sequence_action_invalid_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(invalid_sequence_action='invalid')
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(invalid_sequence_action='Convert')

    def test_invalid_sequence_action_valid_values(self):
        for val in ['ignore', 'fail', 'remove', 'convert', 'convert-ignore', 'convert-remove']:
            self._call_check_inputs(invalid_sequence_action=val)

    def test_alignment_non_bool_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(alignment=1)
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(alignment='true')

    def test_return_list_non_bool_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(return_list=1)

    def test_verbose_non_bool_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(verbose='yes')

    def test_output_filename_non_string_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(output_filename=123)

    def test_output_filename_string_passes(self):
        self._call_check_inputs(output_filename='output.fasta')

    def test_correction_dictionary_non_dict_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(correction_dictionary='not_a_dict')
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(correction_dictionary=[('B', 'N')])

    def test_correction_dictionary_valid_passes(self):
        self._call_check_inputs(correction_dictionary={'B': 'N'})

    def test_ignore_with_expect_unique_incompatible(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(
                duplicate_record_action='ignore',
                expect_unique_header=True,
            )

    def test_ignore_with_expect_unique_false_ok(self):
        self._call_check_inputs(
            duplicate_record_action='ignore',
            expect_unique_header=False,
        )


# ---------------------------------------------------------------------------
# TestInternalParseFasta
# ---------------------------------------------------------------------------
class TestInternalParseFasta:
    """Tests for io.internal_parse_fasta_file and _parse_fasta_all."""

    def test_parse_simple_file(self):
        result = _io.internal_parse_fasta_file(SIMPLE_FILE)
        assert len(result) == 9
        assert isinstance(result, list)
        assert isinstance(result[0], list)
        assert len(result[0]) == 2

    def test_parse_returns_correct_header_and_sequence(self):
        result = _io.internal_parse_fasta_file(SIMPLE_FILE)
        assert result[0][0] == WASL_HEADER
        assert result[0][1] == WASL_SEQ

    def test_sequences_uppercased(self):
        result = _io.internal_parse_fasta_file(SIMPLE_FILE)
        for entry in result:
            assert entry[1] == entry[1].upper()

    def test_file_not_found_raises(self):
        with pytest.raises(ProtfastaException, match='Unable to find file'):
            _io.internal_parse_fasta_file('/nonexistent/path/file.fasta')

    def test_duplicate_headers_with_expect_unique(self):
        with pytest.raises(ProtfastaException, match='duplicate header'):
            _io.internal_parse_fasta_file(DUPLICATE_RECORD_FILE, expect_unique_header=True)

    def test_duplicate_headers_without_expect_unique(self):
        result = _io.internal_parse_fasta_file(DUPLICATE_RECORD_FILE, expect_unique_header=False)
        assert len(result) == 3

    def test_header_parser_applied(self):
        def parser(h):
            return h.split('|')[1] if '|' in h else h
        result = _io.internal_parse_fasta_file(SIMPLE_FILE, header_parser=parser)
        assert result[0][0] == 'O00401'

    def test_parse_fasta_all_empty_content(self):
        result = _io._parse_fasta_all([])
        assert result == []

    def test_parse_fasta_all_blank_lines_skipped(self):
        content = ['>header1\n', 'ACDEF\n', '\n', '>header2\n', 'GHIKL\n']
        result = _io._parse_fasta_all(content)
        assert len(result) == 2
        assert result[0] == ['header1', 'ACDEF']
        assert result[1] == ['header2', 'GHIKL']

    def test_parse_fasta_all_multiline_sequence(self):
        content = ['>header\n', 'ACD\n', 'EFG\n', 'HIK\n']
        result = _io._parse_fasta_all(content)
        assert len(result) == 1
        assert result[0][1] == 'ACDEFGHIK'

    def test_parse_fasta_all_single_entry(self):
        content = ['>header\n', 'ACDEF\n']
        result = _io._parse_fasta_all(content)
        assert len(result) == 1
        assert result[0] == ['header', 'ACDEF']

    def test_verbose_output(self, capsys):
        _io.internal_parse_fasta_file(SIMPLE_FILE, verbose=True)
        captured = capsys.readouterr()
        assert 'Read in file' in captured.out
        assert 'Parsed file' in captured.out


# ---------------------------------------------------------------------------
# TestReadFastaBasic
# ---------------------------------------------------------------------------
class TestReadFastaBasic:
    """Basic read_fasta functionality tests."""

    def test_returns_dict_by_default(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        assert isinstance(result, dict)

    def test_reads_correct_count(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        assert len(result) == 9

    def test_reads_correct_sequence(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        assert result[WASL_HEADER] == WASL_SEQ

    def test_all_values_are_strings(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        for header, seq in result.items():
            assert isinstance(header, str)
            assert isinstance(seq, str)

    def test_all_sequences_nonempty(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        for seq in result.values():
            assert len(seq) > 0

    def test_file_not_found(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta('/nonexistent/file.fasta')


# ---------------------------------------------------------------------------
# TestReadFastaExpectUniqueHeader
# ---------------------------------------------------------------------------
class TestReadFastaExpectUniqueHeader:
    """Tests for the expect_unique_header parameter."""

    def test_true_with_unique_headers(self):
        result = protfasta.read_fasta(SIMPLE_FILE, expect_unique_header=True)
        assert len(result) == 9

    def test_false_with_unique_headers(self):
        result = protfasta.read_fasta(SIMPLE_FILE, expect_unique_header=False)
        assert len(result) == 9

    def test_true_with_duplicate_headers_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(DUPLICATE_RECORD_FILE, expect_unique_header=True)

    def test_non_bool_string_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, expect_unique_header='dog')

    def test_non_bool_int_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, expect_unique_header=1)

    def test_non_bool_none_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, expect_unique_header=None)


# ---------------------------------------------------------------------------
# TestReadFastaHeaderParser
# ---------------------------------------------------------------------------
class TestReadFastaHeaderParser:
    """Tests for the header_parser parameter."""

    def test_truncate_parser(self):
        def truncate(s):
            return s[0:10]
        result = protfasta.read_fasta(SIMPLE_FILE, header_parser=truncate)
        assert len(result) == 9
        assert WASL_HEADER[0:10] in result

    def test_parser_causing_all_same_header_collapses_dict(self):
        def constant(s):
            return "same"
        result = protfasta.read_fasta(
            SIMPLE_FILE,
            header_parser=constant,
            duplicate_sequence_action='ignore',
            expect_unique_header=False,
        )
        assert len(result) == 1  # dict overwrites

    def test_parser_causing_all_same_header_preserved_in_list(self):
        def constant(s):
            return "same"
        result = protfasta.read_fasta(
            SIMPLE_FILE,
            header_parser=constant,
            duplicate_sequence_action='ignore',
            expect_unique_header=False,
            return_list=True,
        )
        assert len(result) == 9

    def test_duplicate_parsed_headers_raise_with_expect_unique(self):
        def constant(s):
            return "same"
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, header_parser=constant)

    def test_no_arg_function_fails(self):
        def no_arg():
            return "hello"
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, header_parser=no_arg)

    def test_function_that_fails_on_test_string(self):
        def bad_parser(s):
            return s.split('|')[1]
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, header_parser=bad_parser, check_header_parser=True)

    def test_check_header_parser_disabled_allows_bad_parser(self):
        def bad_parser(s):
            return s.split('|')[1]
        # Should not raise during input checking
        result = protfasta.read_fasta(SIMPLE_FILE, header_parser=bad_parser, check_header_parser=False)
        assert len(result) == 9

    def test_uniprot_id_parser(self):
        def uniprot_id(s):
            return s.split('|')[1] if '|' in s else s
        result = protfasta.read_fasta(SIMPLE_FILE, header_parser=uniprot_id)
        assert 'O00401' in result
        assert result['O00401'] == WASL_SEQ


# ---------------------------------------------------------------------------
# TestReadFastaDuplicateRecords
# ---------------------------------------------------------------------------
class TestReadFastaDuplicateRecords:
    """Tests for the duplicate_record_action parameter."""

    def test_fail_on_unique_data(self):
        result = protfasta.read_fasta(SIMPLE_FILE, duplicate_record_action='fail')
        assert len(result) == 9

    def test_fail_on_duplicate_data(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(DUPLICATE_RECORD_FILE, duplicate_record_action='fail')

    def test_ignore_requires_expect_unique_false(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(DUPLICATE_RECORD_FILE, duplicate_record_action='ignore')

    def test_ignore_with_expect_unique_false_as_list(self):
        result = protfasta.read_fasta(
            DUPLICATE_RECORD_FILE,
            duplicate_record_action='ignore',
            expect_unique_header=False,
            return_list=True,
        )
        assert len(result) == 3

    def test_ignore_with_expect_unique_false_as_dict(self):
        result = protfasta.read_fasta(
            DUPLICATE_RECORD_FILE,
            duplicate_record_action='ignore',
            expect_unique_header=False,
        )
        # Dict overwrites duplicate key -> 2 unique headers
        assert len(result) == 2

    def test_remove_with_expect_unique_true_raises(self):
        # Still raises because expect_unique_header=True triggers during parsing
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(DUPLICATE_RECORD_FILE, duplicate_record_action='remove')

    def test_remove_with_expect_unique_false(self):
        result = protfasta.read_fasta(
            DUPLICATE_RECORD_FILE,
            duplicate_record_action='remove',
            expect_unique_header=False,
        )
        assert len(result) == 2

    def test_remove_as_list(self):
        result = protfasta.read_fasta(
            DUPLICATE_RECORD_FILE,
            duplicate_record_action='remove',
            expect_unique_header=False,
            return_list=True,
        )
        assert len(result) == 2

    def test_invalid_action_string_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, duplicate_record_action='invalid')


# ---------------------------------------------------------------------------
# TestReadFastaDuplicateSequences
# ---------------------------------------------------------------------------
class TestReadFastaDuplicateSequences:
    """Tests for the duplicate_sequence_action parameter."""

    def test_fail_on_unique_data(self):
        result = protfasta.read_fasta(SIMPLE_FILE, duplicate_sequence_action='fail')
        assert len(result) == 9

    def test_fail_on_duplicate_sequences(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(DUPLICATE_SEQ_FILE, duplicate_sequence_action='fail')

    def test_ignore_keeps_all(self):
        result = protfasta.read_fasta(DUPLICATE_SEQ_FILE, duplicate_sequence_action='ignore')
        assert len(result) == 3

    def test_remove_keeps_first(self):
        result = protfasta.read_fasta(
            DUPLICATE_SEQ_FILE,
            duplicate_sequence_action='remove',
        )
        assert len(result) == 2

    def test_duplicate_seqs_not_duplicate_records(self):
        # Records are NOT duplicates (different headers), so duplicate_record_action='remove'
        # should not remove anything
        result = protfasta.read_fasta(DUPLICATE_SEQ_FILE, duplicate_record_action='remove')
        assert len(result) == 3

    def test_invalid_action_string_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, duplicate_sequence_action='invalid')


# ---------------------------------------------------------------------------
# TestReadFastaInvalidSequences
# ---------------------------------------------------------------------------
class TestReadFastaInvalidSequences:
    """Tests for the invalid_sequence_action parameter."""

    # --- fail ---
    def test_fail_default_on_bad_chars(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(BADCHAR_FILE)

    def test_fail_explicit_on_bad_chars(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='fail')

    def test_fail_on_nonstandard_chars(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action='fail')

    def test_fail_not_triggered_on_valid_data(self):
        result = protfasta.read_fasta(SIMPLE_FILE, invalid_sequence_action='fail')
        assert len(result) == 9

    # --- ignore ---
    def test_ignore_bad_chars(self):
        result = protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='ignore')
        assert len(result) == 4

    def test_ignore_nonstandard_chars(self):
        result = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action='ignore')
        assert len(result) == 4

    # --- remove ---
    def test_remove_all_bad_chars(self):
        result = protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='remove')
        assert len(result) == 0

    def test_remove_all_nonstandard(self):
        result = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action='remove')
        assert len(result) == 0

    # --- convert ---
    def test_convert_nonstandard_succeeds(self):
        result = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action='convert')
        assert len(result) == 4

    def test_convert_bad_chars_fails(self):
        # Bad chars (like '.') aren't in the standard conversion table
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='convert')

    # --- convert-ignore ---
    def test_convert_ignore_nonstandard(self):
        result = protfasta.read_fasta(NONSTANDARD_FILE, invalid_sequence_action='convert-ignore')
        assert len(result) == 4

    def test_convert_ignore_bad_chars(self):
        result = protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='convert-ignore')
        assert len(result) == 4

    # --- convert-remove ---
    def test_convert_remove_all_bad_chars(self):
        result = protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='convert-remove')
        assert len(result) == 0

    def test_convert_remove_fixable(self):
        result = protfasta.read_fasta(FIXABLE_INVALID_FILE, invalid_sequence_action='convert-remove')
        assert len(result) == 1

    # --- invalid action string ---
    def test_invalid_action_string_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(SIMPLE_FILE, invalid_sequence_action='invalid')


# ---------------------------------------------------------------------------
# TestReadFastaAlignment
# ---------------------------------------------------------------------------
class TestReadFastaAlignment:
    """Tests for the alignment parameter."""

    def test_alignment_preserves_dashes(self):
        result = protfasta.read_fasta(ALIGNED_VALID_FILE, alignment=True)
        assert result['Seq1'] == 'A-----CDEFGHIKLMNPQRSTVWY'
        assert result['Seq2'] == 'ACDEFGHIKL-----MNPQRSTVWY'
        assert result['Seq3'] == 'ACDEFGHIKLMNPQRSTVWY-----'

    def test_alignment_reads_all_three(self):
        result = protfasta.read_fasta(ALIGNED_VALID_FILE, alignment=True)
        assert len(result) == 3

    def test_without_alignment_dashes_fail(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(ALIGNED_VALID_FILE)

    def test_non_bool_alignment_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(ALIGNED_VALID_FILE, alignment=1)

    def test_aligned_convertable_fails_by_default(self):
        # Has convertable chars (*, B, Z) that should trigger fail
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(ALIGNED_CONVERTABLE_FILE, alignment=True)

    def test_aligned_convertable_with_convert(self):
        result = protfasta.read_fasta(
            ALIGNED_CONVERTABLE_FILE,
            alignment=True,
            invalid_sequence_action='convert',
        )
        assert result['Seq1'] == 'A-----CDEFGHIKLMNPQRSTVWY'
        assert len(result) == 3

    def test_aligned_unconvertable_fails_with_convert(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(
                ALIGNED_UNCONVERTABLE_FILE,
                alignment=True,
                invalid_sequence_action='convert',
            )

    def test_aligned_unconvertable_fails_by_default(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(ALIGNED_UNCONVERTABLE_FILE, alignment=True)

    def test_aligned_unconvertable_convert_ignore(self):
        result = protfasta.read_fasta(
            ALIGNED_UNCONVERTABLE_FILE,
            alignment=True,
            invalid_sequence_action='convert-ignore',
        )
        assert result['Seq2'] == 'ACDEFGHIKL-----MNPQRSTVWYN'

    def test_aligned_unconvertable_remove_all(self):
        result = protfasta.read_fasta(
            ALIGNED_UNCONVERTABLE_FILE,
            alignment=True,
            invalid_sequence_action='remove',
        )
        assert len(result) == 0

    def test_aligned_valid_remove_keeps_all(self):
        result = protfasta.read_fasta(
            ALIGNED_VALID_FILE,
            alignment=True,
            invalid_sequence_action='remove',
        )
        assert len(result) == 3

    def test_without_alignment_remove_removes_dashed(self):
        # Without alignment flag, dashes are invalid -> all removed
        result = protfasta.read_fasta(
            ALIGNED_VALID_FILE,
            invalid_sequence_action='remove',
        )
        assert len(result) == 0


# ---------------------------------------------------------------------------
# TestReadFastaReturnList
# ---------------------------------------------------------------------------
class TestReadFastaReturnList:
    """Tests for the return_list parameter."""

    def test_returns_dict_by_default(self):
        result = protfasta.read_fasta(SIMPLE_FILE)
        assert isinstance(result, dict)

    def test_returns_list_when_true(self):
        result = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        assert isinstance(result, list)

    def test_list_elements_are_pairs(self):
        result = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        for entry in result:
            assert len(entry) == 2

    def test_list_preserves_order(self):
        result = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        assert result[0][0] == WASL_HEADER
        assert result[0][1] == WASL_SEQ

    def test_list_preserves_duplicate_records(self):
        result = protfasta.read_fasta(
            DUPLICATE_RECORD_FILE,
            duplicate_record_action='ignore',
            return_list=True,
            expect_unique_header=False,
        )
        assert len(result) == 3

    def test_list_count_matches_dict_count(self):
        d = protfasta.read_fasta(SIMPLE_FILE)
        l = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        assert len(d) == len(l)


# ---------------------------------------------------------------------------
# TestReadFastaOutputFile
# ---------------------------------------------------------------------------
class TestReadFastaOutputFile:
    """Tests for the output_filename parameter."""

    def test_output_file_written(self, tmp_path):
        outfile = str(tmp_path / 'output.fasta')
        result = protfasta.read_fasta(SIMPLE_FILE, output_filename=outfile)
        assert os.path.exists(outfile)
        # Read back and verify
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 9
        for k in result:
            assert readback[k] == result[k]

    def test_output_file_after_filtering(self, tmp_path):
        outfile = str(tmp_path / 'filtered.fasta')
        result = protfasta.read_fasta(
            DUPLICATE_SEQ_FILE,
            duplicate_sequence_action='remove',
            output_filename=outfile,
        )
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 2


# ---------------------------------------------------------------------------
# TestReadFastaCorrectionDictionary
# ---------------------------------------------------------------------------
class TestReadFastaCorrectionDictionary:
    """Tests for the correction_dictionary parameter."""

    def test_correction_dict_without_convert_fails(self):
        # Passing a correction dictionary without requesting conversion should fail
        # because the file has invalid sequences and default action is 'fail'
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(NONSTANDARD_FILE, correction_dictionary={'.': 'A'})

    def test_custom_dict_overrides_default(self):
        # Standard dict maps X->G, but custom overrides all conversions
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(
                NONSTANDARD_FILE,
                correction_dictionary={'.': 'A'},
                invalid_sequence_action='convert',
            )

    def test_custom_dict_for_bad_chars(self):
        cd = {'.': 'A', '-': 'C'}
        result = protfasta.read_fasta(
            BADCHAR_FILE,
            correction_dictionary=cd,
            invalid_sequence_action='convert',
        )
        assert len(result) == 4

    def test_incomplete_custom_dict_fails(self):
        # Only maps '.', but file also has '-'
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(
                BADCHAR_FILE,
                correction_dictionary={'.': 'A'},
                invalid_sequence_action='convert',
            )

    def test_custom_dict_with_convert_ignore(self):
        cd = {'.': 'A'}
        result = protfasta.read_fasta(
            BADCHAR_FILE,
            correction_dictionary=cd,
            invalid_sequence_action='convert-ignore',
        )
        assert len(result) == 4


# ---------------------------------------------------------------------------
# TestWriteFasta
# ---------------------------------------------------------------------------
class TestWriteFasta:
    """Tests for write_fasta function."""

    def test_write_dict_roundtrip(self, tmp_path):
        original = protfasta.read_fasta(SIMPLE_FILE)
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta(original, outfile)
        readback = protfasta.read_fasta(outfile)
        for k in original:
            assert readback[k] == original[k]

    def test_write_list_roundtrip(self, tmp_path):
        original = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta(original, outfile)
        readback = protfasta.read_fasta(outfile, return_list=True)
        for idx in range(len(original)):
            assert readback[idx][0] == original[idx][0]
            assert readback[idx][1] == original[idx][1]

    def test_write_and_read_preserves_count(self, tmp_path):
        original = protfasta.read_fasta(SIMPLE_FILE)
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta(original, outfile)
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == len(original)

    def test_append_mode(self, tmp_path):
        original = protfasta.read_fasta(SIMPLE_FILE)
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta(original, outfile)

        added = {'added_sequence': 'ASPAPSPAPSPAPSPAS'}
        protfasta.write_fasta(added, outfile, append_to_fasta=True)

        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 10
        assert readback['added_sequence'] == 'ASPAPSPAPSPAPSPAS'
        for k in original:
            assert readback[k] == original[k]

    def test_append_to_nonexistent_creates(self, tmp_path):
        outfile = str(tmp_path / 'new.fasta')
        data = {'header1': 'ACDEF'}
        protfasta.write_fasta(data, outfile, append_to_fasta=True)
        readback = protfasta.read_fasta(outfile)
        assert readback['header1'] == 'ACDEF'

    def test_overwrite_mode(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta({'h1': 'ACDEF'}, outfile)
        protfasta.write_fasta({'h2': 'GHIKL'}, outfile)  # overwrites
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 1
        assert 'h2' in readback

    def test_linelength_default(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 120
        protfasta.write_fasta({'header': seq}, outfile, linelength=60)
        with open(outfile) as f:
            lines = f.readlines()
        # Header line + 2 sequence lines (60+60) + possible trailing newline
        seq_lines = [l for l in lines if not l.startswith('>') and l.strip()]
        assert len(seq_lines) == 2
        assert len(seq_lines[0].strip()) == 60

    def test_linelength_none_no_wrap(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 200
        protfasta.write_fasta({'header': seq}, outfile, linelength=None)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [l for l in lines if not l.startswith('>') and l.strip()]
        assert len(seq_lines) == 1
        assert len(seq_lines[0].strip()) == 200

    def test_linelength_false_no_wrap(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 200
        protfasta.write_fasta({'header': seq}, outfile, linelength=False)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [l for l in lines if not l.startswith('>') and l.strip()]
        assert len(seq_lines) == 1

    def test_linelength_very_short_clamped_to_5(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 20
        protfasta.write_fasta({'header': seq}, outfile, linelength=2)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [l for l in lines if not l.startswith('>') and l.strip()]
        # With linelength=5: 20/5 = 4 lines
        assert len(seq_lines) == 4

    def test_write_empty_sequence_raises(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta({'header': ''}, outfile)

    def test_write_list_bad_element_raises(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta([['header_only']], outfile)

    def test_write_multiple_sequences(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        data = {'h1': 'ACDEF', 'h2': 'GHIKL', 'h3': 'MNPQR'}
        protfasta.write_fasta(data, outfile)
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 3
        for k in data:
            assert readback[k] == data[k]


# ---------------------------------------------------------------------------
# TestEndToEndCombinations
# ---------------------------------------------------------------------------
class TestEndToEndCombinations:
    """Tests for combinations of parameters and edge cases."""

    def test_verbose_flag(self, capsys):
        protfasta.read_fasta(
            DUPLICATE_SEQ_FILE,
            duplicate_sequence_action='remove',
            verbose=True,
        )
        captured = capsys.readouterr()
        assert 'duplicate sequences' in captured.out.lower() or 'Removed' in captured.out

    def test_convert_remove_then_read_back(self, tmp_path):
        outfile = str(tmp_path / 'output.fasta')
        result = protfasta.read_fasta(
            FIXABLE_INVALID_FILE,
            invalid_sequence_action='convert-remove',
            output_filename=outfile,
        )
        assert len(result) == 1
        # Read back - should be clean standard amino acids
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 1

    def test_alignment_write_roundtrip(self, tmp_path):
        outfile = str(tmp_path / 'aligned.fasta')
        result = protfasta.read_fasta(ALIGNED_VALID_FILE, alignment=True)
        protfasta.write_fasta(result, outfile)
        readback = protfasta.read_fasta(outfile, alignment=True)
        for k in result:
            assert readback[k] == result[k]

    def test_header_parser_with_duplicate_removal(self):
        def first_word(s):
            return s.split()[0]
        result = protfasta.read_fasta(
            SIMPLE_FILE,
            header_parser=first_word,
            duplicate_sequence_action='remove',
        )
        assert len(result) == 9

    def test_return_list_with_duplicate_removal(self):
        result = protfasta.read_fasta(
            DUPLICATE_SEQ_FILE,
            duplicate_sequence_action='remove',
            return_list=True,
        )
        assert isinstance(result, list)
        assert len(result) == 2

    def test_read_write_read_consistency(self, tmp_path):
        """Read -> write -> read should give identical data."""
        outfile = str(tmp_path / 'roundtrip.fasta')
        original = protfasta.read_fasta(SIMPLE_FILE)
        protfasta.write_fasta(original, outfile)
        readback = protfasta.read_fasta(outfile)
        assert original == readback

    def test_read_write_read_list_consistency(self, tmp_path):
        """Read (list) -> write -> read (list) should give identical data."""
        outfile = str(tmp_path / 'roundtrip.fasta')
        original = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        protfasta.write_fasta(original, outfile)
        readback = protfasta.read_fasta(outfile, return_list=True)
        assert len(original) == len(readback)
        for o, r in zip(original, readback):
            assert o[0] == r[0]
            assert o[1] == r[1]

    def test_large_linelength_roundtrip(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        original = protfasta.read_fasta(SIMPLE_FILE)
        protfasta.write_fasta(original, outfile, linelength=200)
        readback = protfasta.read_fasta(outfile)
        assert original == readback

    def test_multiple_appends(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta({'h1': 'ACDEF'}, outfile)
        protfasta.write_fasta({'h2': 'GHIKL'}, outfile, append_to_fasta=True)
        protfasta.write_fasta({'h3': 'MNPQR'}, outfile, append_to_fasta=True)
        readback = protfasta.read_fasta(outfile)
        assert len(readback) == 3
        assert readback['h1'] == 'ACDEF'
        assert readback['h2'] == 'GHIKL'
        assert readback['h3'] == 'MNPQR'
