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
- TestNonASCIISequences: Non-ASCII characters in sequence data
- TestPathLikeArguments: pathlib.Path accepted wherever a filename is
- TestFileOpenErrors: OS-level file errors surface as ProtfastaException
- TestEndToEndCombinations: Combined parameter interactions
- TestReadFastaStream: read_fasta_stream streaming parser
- TestDuplicateRecordsMultiSequenceHeader: duplicate records under a header
  that also carries other sequences
- TestEncoding: byte-order marks and non-UTF-8 bytes
- TestParserEdgeCases: line-level edge cases in the parsing engine
- TestConverter: the shared sequence converter
- TestWriteFastaAtomicity: write_fasta never leaves a partial file
- TestStreamRobustness: file-handling edge cases when streaming
- TestReadFastaStreamParity: read_fasta and read_fasta_stream agree
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

    def test_first_invalid_char_reported_in_sequence_order(self):
        valid, info = _utilities.check_sequence_is_valid('ACD*E-F.G')
        assert valid is False
        assert info == '*'

    def test_non_ascii_char_invalid(self):
        # exercises the non-ASCII (str.translate) path
        valid, info = _utilities.check_sequence_is_valid('ACDÉF')
        assert valid is False
        assert info == 'É'

    def test_non_ascii_with_alignment(self):
        assert _utilities.check_sequence_is_valid('A--CDÉ', alignment=True) == (False, 'É')
        assert _utilities.check_sequence_is_valid('A--CD\u00e9'.upper(), alignment=True)[0] is False

    def test_lone_surrogate_invalid(self):
        # what an undecodable byte looks like after surrogateescape decoding
        valid, info = _utilities.check_sequence_is_valid('ACD\udce9')
        assert valid is False
        assert info == '\udce9'

    def test_all_valid_ascii_and_non_ascii_paths_agree(self):
        # the ASCII fast path and the str fallback must give the same answer
        for seq in ['ACDEFGHIKLMNPQRSTVWY', 'ACD-', 'ACDB', '', 'X']:
            expected = _utilities.check_sequence_is_valid(seq)
            # force the non-ASCII path by appending and stripping a marker
            valid, info = _utilities.check_sequence_is_valid(seq + 'É')
            assert valid is False
            if expected[0]:
                assert info == 'É'
            else:
                assert info == expected[1]


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

    def test_multi_character_key(self):
        cd = {'XX': 'G', 'B': 'N'}
        assert _utilities.convert_to_valid('AXXB', correction_dictionary=cd) == 'AGN'

    def test_multi_character_key_applied_in_dictionary_order(self):
        cd = {'AB': 'C', 'CC': 'D'}
        assert _utilities.convert_to_valid('ABC', correction_dictionary=cd) == 'D'

    def test_empty_key_raises(self):
        # an empty key would make str.replace insert between every residue
        with pytest.raises(ProtfastaException):
            _utilities.convert_to_valid('ACD', correction_dictionary={'': 'X'})

    def test_non_string_value_raises(self):
        with pytest.raises(ProtfastaException):
            _utilities.convert_to_valid('ACD', correction_dictionary={'D': 5})

    def test_non_string_key_raises(self):
        with pytest.raises(ProtfastaException):
            _utilities.convert_to_valid('ACD', correction_dictionary={5: 'D'})

    def test_clean_sequence_returned_unchanged(self):
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        assert _utilities.convert_to_valid(seq) is seq
        assert _utilities.convert_to_valid(seq, alignment=True) is seq
        assert _utilities.convert_to_valid(seq, correction_dictionary={'.': 'A'}) is seq

    def test_non_ascii_untouched_by_builtin_table(self):
        assert _utilities.convert_to_valid('ACDÉ') == 'ACDÉ'


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

    def test_fail_on_duplicates_after_header_reused_with_other_sequence(self):
        # Regression: the old lookup only remembered the FIRST sequence seen
        # for a header, so once 'h1' appeared with a second sequence, a
        # repeat of that second sequence slipped through undetected.
        data = [['h1', 'ACDEF'], ['h1', 'GHIKL'], ['h1', 'GHIKL']]
        with pytest.raises(ProtfastaException):
            _utilities.fail_on_duplicates(data)

    def test_remove_duplicates_after_header_reused_with_other_sequence(self):
        data = [['h1', 'ACDEF'], ['h1', 'GHIKL'], ['h1', 'GHIKL'], ['h1', 'ACDEF']]
        result = _utilities.remove_duplicates(data)
        assert result == [['h1', 'ACDEF'], ['h1', 'GHIKL']]

    def test_record_hash_distinguishes_header_sequence_boundary(self):
        # The header/sequence split must be part of the key
        assert _utilities._record_hash('AB', 'CD') != _utilities._record_hash('A', 'BCD')
        assert _utilities._record_hash('h', 'ACD') == _utilities._record_hash('h', 'ACD')

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

    def test_correction_dictionary_empty_key_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(correction_dictionary={'': 'N'})

    def test_correction_dictionary_non_string_value_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(correction_dictionary={'B': None})

    def test_correction_dictionary_multi_char_key_passes(self):
        self._call_check_inputs(correction_dictionary={'XX': 'G'})

    def test_header_parser_non_callable_fails_even_when_check_disabled(self):
        # check_header_parser only controls the smoke test; a non-callable
        # can never work and used to surface as a TypeError mid-parse
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(header_parser='not_a_function', check_header_parser=False)

    def test_check_header_parser_non_bool_fails(self):
        with pytest.raises(ProtfastaException):
            self._call_check_inputs(check_header_parser='yes')

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

    def test_parser_raising_mid_parse_is_a_protfasta_exception(self, tmp_path):
        # A parser that passes the smoke test but blows up on a real header
        # must still surface as a ProtfastaException, naming the header.
        f = tmp_path / 'in.fasta'
        f.write_text('>sp|P1|OK\nACDE\n>no_pipe_here\nFGHI\n')

        def accession(s):
            return s.split('|')[1] if '|' not in 'this test string should work' or True else s

        def strict(s):
            return s.split('|')[1]

        with pytest.raises(ProtfastaException, match='no_pipe_here'):
            protfasta.read_fasta(str(f), header_parser=strict, check_header_parser=False)
        with pytest.raises(ProtfastaException, match='no_pipe_here'):
            list(protfasta.read_fasta_stream(str(f), header_parser=strict, check_header_parser=False))
        del accession

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
        as_dict = protfasta.read_fasta(SIMPLE_FILE)
        as_list = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        assert len(as_dict) == len(as_list)


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
        assert readback == result


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
        seq_lines = [line for line in lines if not line.startswith('>') and line.strip()]
        assert len(seq_lines) == 2
        assert len(seq_lines[0].strip()) == 60

    def test_linelength_none_no_wrap(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 200
        protfasta.write_fasta({'header': seq}, outfile, linelength=None)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [line for line in lines if not line.startswith('>') and line.strip()]
        assert len(seq_lines) == 1
        assert len(seq_lines[0].strip()) == 200

    def test_linelength_false_no_wrap(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 200
        protfasta.write_fasta({'header': seq}, outfile, linelength=False)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [line for line in lines if not line.startswith('>') and line.strip()]
        assert len(seq_lines) == 1

    def test_linelength_very_short_clamped_to_5(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        seq = 'A' * 20
        protfasta.write_fasta({'header': seq}, outfile, linelength=2)
        with open(outfile) as f:
            lines = f.readlines()
        seq_lines = [line for line in lines if not line.startswith('>') and line.strip()]
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

    def test_unsupported_container_raises(self, tmp_path):
        # A tuple of pairs is not supported - this should be a clean
        # ProtfastaException rather than an UnboundLocalError.
        outfile = str(tmp_path / 'test.fasta')
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta((('h1', 'ACDEF'),), outfile)

    def test_generator_input_raises(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta((x for x in [['h1', 'ACDEF']]), outfile)

    def test_numeric_string_linelength_accepted(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta({'header': 'A' * 120}, outfile, linelength='60')
        with open(outfile) as f:
            seq_lines = [line for line in f if not line.startswith('>') and line.strip()]
        assert len(seq_lines) == 2

    def test_non_numeric_linelength_raises(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta({'header': 'ACDEF'}, outfile, linelength='sixty')

    def test_exact_output_format(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta({'h1': 'A' * 12, 'h2': 'CDE'}, outfile, linelength=5)
        with open(outfile) as f:
            assert f.read() == '>h1\nAAAAA\nAAAAA\nAA\n\n>h2\nCDE\n\n'

    def test_exact_output_format_single_line(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta([['h1', 'A' * 12]], outfile, linelength=None)
        with open(outfile) as f:
            assert f.read() == '>h1\nAAAAAAAAAAAA\n\n'

    def test_dict_and_list_input_write_identical_files(self, tmp_path):
        d = protfasta.read_fasta(SIMPLE_FILE)
        as_list = [[k, v] for k, v in d.items()]
        f1 = tmp_path / 'a.fasta'
        f2 = tmp_path / 'b.fasta'
        protfasta.write_fasta(d, f1)
        protfasta.write_fasta(as_list, f2)
        assert f1.read_bytes() == f2.read_bytes()

    def test_linelength_zero_no_wrap(self, tmp_path):
        outfile = str(tmp_path / 'test.fasta')
        protfasta.write_fasta({'header': 'A' * 200}, outfile, linelength=0)
        with open(outfile) as f:
            seq_lines = [line for line in f if not line.startswith('>') and line.strip()]
        assert len(seq_lines) == 1


# ---------------------------------------------------------------------------
# TestNonASCIISequences
# ---------------------------------------------------------------------------
class TestNonASCIISequences:
    """Non-ASCII characters in sequence data must not crash the parser.

    Duplicate detection hashes every sequence *before* invalid-residue
    handling gets a chance to flag or strip it, so the hashing step has to
    tolerate arbitrary characters.
    """

    @pytest.fixture
    def nonascii_file(self, tmp_path):
        f = tmp_path / 'nonascii.fasta'
        f.write_text('>h1\nACDEÉF\n>h2\nGHIKL\n', encoding='utf-8')
        return str(f)

    def test_seq_hash_handles_non_ascii(self):
        digest = _utilities._seq_hash('ACDEÉF')
        assert isinstance(digest, bytes)
        assert len(digest) == 16

    def test_seq_hash_distinguishes_non_ascii(self):
        assert _utilities._seq_hash('ACDÉ') != _utilities._seq_hash('ACDÊ')

    def test_ignore_keeps_non_ascii_record(self, nonascii_file):
        result = protfasta.read_fasta(nonascii_file, invalid_sequence_action='ignore')
        assert len(result) == 2
        assert result['h1'] == 'ACDEÉF'

    def test_fail_raises_protfasta_exception(self, nonascii_file):
        # Must be a ProtfastaException, not a UnicodeEncodeError
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(nonascii_file)

    def test_remove_drops_non_ascii_record(self, nonascii_file):
        result = protfasta.read_fasta(nonascii_file, invalid_sequence_action='remove')
        assert len(result) == 1
        assert 'h2' in result

    def test_duplicate_sequence_remove_with_non_ascii(self, nonascii_file):
        result = protfasta.read_fasta(
            nonascii_file,
            invalid_sequence_action='ignore',
            duplicate_sequence_action='remove',
        )
        assert len(result) == 2

    def test_streaming_handles_non_ascii(self, nonascii_file):
        streamed = list(protfasta.read_fasta_stream(
            nonascii_file,
            invalid_sequence_action='ignore',
            duplicate_sequence_action='remove',
            silence_warnings=True,
        ))
        assert len(streamed) == 2


# ---------------------------------------------------------------------------
# TestPathLikeArguments
# ---------------------------------------------------------------------------
class TestPathLikeArguments:
    """pathlib.Path objects are accepted anywhere a filename is expected."""

    def test_read_fasta_accepts_path(self):
        result = protfasta.read_fasta(Path(SIMPLE_FILE))
        assert len(result) == 9

    def test_read_fasta_output_filename_accepts_path(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        protfasta.read_fasta(SIMPLE_FILE, output_filename=outfile)
        assert outfile.exists()

    def test_write_fasta_accepts_path(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        protfasta.write_fasta({'h1': 'ACDEF'}, outfile)
        assert protfasta.read_fasta(outfile)['h1'] == 'ACDEF'

    def test_read_fasta_stream_accepts_path(self):
        streamed = list(protfasta.read_fasta_stream(Path(SIMPLE_FILE)))
        assert len(streamed) == 9

    def test_stream_output_filename_accepts_path(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        list(protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=outfile))
        assert len(protfasta.read_fasta(outfile)) == 9

    def test_integer_filename_rejected(self):
        # open() would treat an int as a file descriptor and silently read
        # from it; we want a clear error instead.
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(0)

    def test_integer_filename_rejected_streaming(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta_stream(0)


# ---------------------------------------------------------------------------
# TestFileOpenErrors
# ---------------------------------------------------------------------------
class TestFileOpenErrors:
    """Every file-open failure mode surfaces as a ProtfastaException."""

    def test_directory_instead_of_file(self, tmp_path):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta(str(tmp_path))

    def test_directory_instead_of_file_streaming(self, tmp_path):
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream(str(tmp_path)))

    def test_missing_file_message(self):
        with pytest.raises(ProtfastaException, match='Unable to find file'):
            protfasta.read_fasta('/nonexistent/path/file.fasta')


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


# ---------------------------------------------------------------------------
# TestReadFastaStream
# ---------------------------------------------------------------------------
class TestReadFastaStream:
    """Tests for read_fasta_stream, the streaming counterpart to read_fasta."""

    def test_is_callable(self):
        assert callable(protfasta.read_fasta_stream)

    def test_returns_generator(self):
        import types
        result = protfasta.read_fasta_stream(SIMPLE_FILE)
        assert isinstance(result, types.GeneratorType)

    def test_yields_tuples_by_default(self):
        for record in protfasta.read_fasta_stream(SIMPLE_FILE):
            assert isinstance(record, tuple)
            assert len(record) == 2
            break

    def test_yields_lists_with_return_list(self):
        for record in protfasta.read_fasta_stream(SIMPLE_FILE, return_list=True):
            assert isinstance(record, list)
            assert len(record) == 2
            break

    def test_count_matches_read_fasta(self):
        ref = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        streamed = list(protfasta.read_fasta_stream(SIMPLE_FILE))
        assert len(streamed) == len(ref)

    def test_content_parity_with_read_fasta(self):
        ref = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        streamed = list(protfasta.read_fasta_stream(SIMPLE_FILE))
        for (h_ref, s_ref), (h_str, s_str) in zip(ref, streamed):
            assert h_ref == h_str
            assert s_ref == s_str

    def test_preserves_order(self):
        streamed = list(protfasta.read_fasta_stream(SIMPLE_FILE))
        assert streamed[0][0] == WASL_HEADER
        assert streamed[0][1] == WASL_SEQ

    def test_header_parser_applied(self):
        def first_word(s):
            return s.split()[0]
        streamed = list(protfasta.read_fasta_stream(SIMPLE_FILE, header_parser=first_word))
        assert streamed[0][0] == WASL_HEADER.split()[0]

    def test_invalid_action_fail_raises(self):
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream(BADCHAR_FILE))

    def test_invalid_action_remove(self):
        ref = protfasta.read_fasta(BADCHAR_FILE, invalid_sequence_action='remove', return_list=True)
        streamed = list(protfasta.read_fasta_stream(BADCHAR_FILE, invalid_sequence_action='remove'))
        assert len(streamed) == len(ref)

    def test_invalid_action_convert_parity(self):
        ref = protfasta.read_fasta(FIXABLE_INVALID_FILE, invalid_sequence_action='convert', return_list=True)
        streamed = list(protfasta.read_fasta_stream(FIXABLE_INVALID_FILE, invalid_sequence_action='convert'))
        assert [list(r) for r in streamed] == ref

    def test_convert_remove_parity(self):
        ref = protfasta.read_fasta(FIXABLE_INVALID_FILE, invalid_sequence_action='convert-remove', return_list=True)
        streamed = list(protfasta.read_fasta_stream(FIXABLE_INVALID_FILE, invalid_sequence_action='convert-remove'))
        assert [list(r) for r in streamed] == ref

    # Note the tests below opt into memory-growing duplicate checks but are not
    # testing the warning that accompanies them, so they pass
    # silence_warnings=True to keep the test output clean. The warning itself is
    # covered by the dedicated tests at the end of this class.

    def test_duplicate_sequence_remove_parity(self):
        ref = protfasta.read_fasta(DUPLICATE_SEQ_FILE, duplicate_sequence_action='remove', return_list=True)
        streamed = list(protfasta.read_fasta_stream(DUPLICATE_SEQ_FILE,
                                                    duplicate_sequence_action='remove',
                                                    silence_warnings=True))
        assert [list(r) for r in streamed] == ref

    def test_duplicate_sequence_fail_raises(self):
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream(DUPLICATE_SEQ_FILE,
                                             duplicate_sequence_action='fail',
                                             silence_warnings=True))

    def test_duplicate_record_fail_raises(self):
        # duplicate_record_action defaults to 'ignore' (flat memory), so opt in
        # to 'fail' explicitly to exercise duplicate-record detection.
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream(DUPLICATE_RECORD_FILE,
                                             expect_unique_header=False,
                                             duplicate_record_action='fail',
                                             silence_warnings=True))

    def test_eager_validation_bad_kwarg(self):
        # Bad keyword must raise at call time, before any iteration begins.
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta_stream(SIMPLE_FILE, invalid_sequence_action='nonsense')

    def test_missing_file_raises(self):
        # File-open errors surface when iteration begins.
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream('does_not_exist_12345.fasta'))

    def test_output_filename_tee_parity(self, tmp_path):
        outfile = str(tmp_path / 'streamed.fasta')
        # Fully consume the stream so the tee'd file is complete.
        list(protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=outfile))
        ref = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        readback = protfasta.read_fasta(outfile, return_list=True)
        assert readback == ref

    def test_output_filename_same_as_input_raises(self):
        with pytest.raises(ProtfastaException):
            protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=SIMPLE_FILE)

    def test_alignment_parity(self):
        ref = protfasta.read_fasta(ALIGNED_VALID_FILE, alignment=True, return_list=True)
        streamed = list(protfasta.read_fasta_stream(ALIGNED_VALID_FILE, alignment=True))
        assert [list(r) for r in streamed] == ref

    def test_verbose_summary_at_exhaustion(self, capsys):
        stream = protfasta.read_fasta_stream(
            DUPLICATE_SEQ_FILE,
            duplicate_sequence_action='remove',
            verbose=True,
            silence_warnings=True,
        )
        # Nothing summarised until the stream is consumed.
        list(stream)
        captured = capsys.readouterr()
        assert 'duplicate sequences' in captured.out.lower() or 'Streamed' in captured.out

    # -- memory-growth warning / flat-memory mode -----------------------------

    def test_flat_by_default_no_warning(self):
        """The defaults are memory-flat, so no warning is emitted."""
        import warnings
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            protfasta.read_fasta_stream(SIMPLE_FILE)
        assert not [w for w in caught if "O(number of records)" in str(w.message)]

    def test_warns_when_check_enabled(self):
        """Opting into a duplicate/uniqueness check emits a one-time memory warning."""
        with pytest.warns(UserWarning, match=r"O\(number of records\)"):
            protfasta.read_fasta_stream(SIMPLE_FILE, expect_unique_header=True)

    def test_warns_lists_the_enabled_checks(self):
        """The warning names the specific options that cause the growth."""
        with pytest.warns(UserWarning, match="duplicate_sequence_action"):
            protfasta.read_fasta_stream(
                SIMPLE_FILE,
                expect_unique_header=False,
                duplicate_record_action='ignore',
                duplicate_sequence_action='fail',
            )

    def test_no_warning_when_flat(self):
        """With every check disabled, streaming emits no memory warning."""
        import warnings
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            protfasta.read_fasta_stream(
                SIMPLE_FILE,
                expect_unique_header=False,
                duplicate_record_action='ignore',
                duplicate_sequence_action='ignore',
            )
        assert not [w for w in caught if "O(number of records)" in str(w.message)]

    def test_silence_warnings_suppresses(self):
        """silence_warnings=True suppresses the memory warning even with checks on."""
        import warnings
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            protfasta.read_fasta_stream(SIMPLE_FILE, silence_warnings=True)
        assert not [w for w in caught if "O(number of records)" in str(w.message)]

    def test_flat_options_content_parity(self):
        """Flat streaming yields the same records as read_fasta."""
        import warnings
        ref = protfasta.read_fasta(SIMPLE_FILE, return_list=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            streamed = list(protfasta.read_fasta_stream(
                SIMPLE_FILE,
                expect_unique_header=False,
                duplicate_record_action='ignore',
                duplicate_sequence_action='ignore',
            ))
        assert [list(r) for r in streamed] == ref


# ---------------------------------------------------------------------------
# TestDuplicateRecordsMultiSequenceHeader
# ---------------------------------------------------------------------------
class TestDuplicateRecordsMultiSequenceHeader:
    """A header that appears with several different sequences, one of which
    is then repeated.  The repeat is a duplicate *record* and must be caught
    by read_fasta exactly as read_fasta_stream already caught it."""

    @pytest.fixture
    def multi_seq_file(self, tmp_path):
        f = tmp_path / 'multi.fasta'
        f.write_text('>h\nACDE\n>h\nFGHI\n>h\nFGHI\n>h\nACDE\n>other\nKLMN\n')
        return str(f)

    def test_read_fasta_fail_detects_repeat(self, multi_seq_file):
        with pytest.raises(ProtfastaException, match='duplicate entries'):
            protfasta.read_fasta(multi_seq_file, expect_unique_header=False, duplicate_record_action='fail')

    def test_read_fasta_remove_keeps_first_of_each(self, multi_seq_file):
        result = protfasta.read_fasta(
            multi_seq_file,
            expect_unique_header=False,
            duplicate_record_action='remove',
            return_list=True,
        )
        assert result == [['h', 'ACDE'], ['h', 'FGHI'], ['other', 'KLMN']]

    def test_stream_agrees_with_read_fasta(self, multi_seq_file):
        ref = protfasta.read_fasta(
            multi_seq_file,
            expect_unique_header=False,
            duplicate_record_action='remove',
            return_list=True,
        )
        streamed = list(protfasta.read_fasta_stream(
            multi_seq_file,
            duplicate_record_action='remove',
            return_list=True,
            silence_warnings=True,
        ))
        assert streamed == ref

    def test_stream_fail_detects_repeat(self, multi_seq_file):
        with pytest.raises(ProtfastaException, match='duplicate entries'):
            list(protfasta.read_fasta_stream(multi_seq_file, duplicate_record_action='fail', silence_warnings=True))

    def test_unique_headers_make_record_check_a_no_op(self):
        # With expect_unique_header=True the record pass is skipped as
        # redundant; the result must be identical to running it.
        a = protfasta.read_fasta(SIMPLE_FILE, duplicate_record_action='remove', return_list=True)
        b = protfasta.read_fasta(SIMPLE_FILE, expect_unique_header=False, duplicate_record_action='remove', return_list=True)
        assert a == b

    def test_unique_headers_verbose_still_reports_record_stage(self, capsys):
        protfasta.read_fasta(SIMPLE_FILE, duplicate_record_action='remove', verbose=True)
        assert 'duplicate records' in capsys.readouterr().out

    def test_stream_warning_omits_record_action_when_headers_unique(self):
        # the record check is skipped when headers are unique, so it must
        # not be blamed for memory growth in the warning
        with pytest.warns(UserWarning) as caught:
            protfasta.read_fasta_stream(SIMPLE_FILE, expect_unique_header=True, duplicate_record_action='remove')
        msgs = ' '.join(str(w.message) for w in caught)
        assert 'expect_unique_header=True' in msgs
        # the remedy text always mentions duplicate_record_action='ignore';
        # what must be absent is the *enabled* action being blamed
        assert "duplicate_record_action='remove'" not in msgs


# ---------------------------------------------------------------------------
# TestEncoding
# ---------------------------------------------------------------------------
class TestEncoding:
    """Byte-order marks and non-UTF-8 bytes must never lose data or escape
    as a UnicodeError."""

    @pytest.fixture
    def bom_file(self, tmp_path):
        f = tmp_path / 'bom.fasta'
        f.write_bytes(b'\xef\xbb\xbf>h1\nACDE\n>h2\nFGHI\n')
        return str(f)

    @pytest.fixture
    def latin1_file(self, tmp_path):
        # 'prot\xe9ine' is Latin-1, not valid UTF-8
        f = tmp_path / 'latin1.fasta'
        f.write_bytes(b'>prot\xe9ine A\nACDE\n>plain\nFGHI\n')
        return str(f)

    @pytest.fixture
    def bad_byte_in_seq_file(self, tmp_path):
        f = tmp_path / 'badseq.fasta'
        f.write_bytes(b'>h1\nAC\xe9DE\n>h2\nFGHI\n')
        return str(f)

    def test_bom_first_record_kept(self, bom_file):
        # Regression: the BOM used to hide the first '>' and silently drop h1
        assert protfasta.read_fasta(bom_file) == {'h1': 'ACDE', 'h2': 'FGHI'}

    def test_bom_first_record_kept_streaming(self, bom_file):
        assert list(protfasta.read_fasta_stream(bom_file)) == [('h1', 'ACDE'), ('h2', 'FGHI')]

    def test_bom_not_stripped_from_later_records(self, tmp_path):
        # only a *leading* BOM is a BOM; elsewhere U+FEFF is just data
        f = tmp_path / 'inner.fasta'
        f.write_bytes(b'>h1\nACDE\n>\xef\xbb\xbfh2\nFGHI\n')
        result = protfasta.read_fasta(str(f))
        assert '\ufeffh2' in result

    def test_latin1_header_reads_without_decode_error(self, latin1_file):
        result = protfasta.read_fasta(latin1_file)
        assert len(result) == 2
        assert 'plain' in result

    def test_latin1_header_roundtrips_byte_for_byte(self, latin1_file, tmp_path):
        result = protfasta.read_fasta(latin1_file, return_list=True)
        out = tmp_path / 'out.fasta'
        protfasta.write_fasta(result, out)
        assert b'>prot\xe9ine A\n' in out.read_bytes()
        assert protfasta.read_fasta(str(out), return_list=True) == result

    def test_latin1_header_streams_and_tees(self, latin1_file, tmp_path):
        out = tmp_path / 'out.fasta'
        streamed = list(protfasta.read_fasta_stream(latin1_file, output_filename=out))
        assert len(streamed) == 2
        assert b'>prot\xe9ine A\n' in out.read_bytes()

    def test_bad_byte_in_sequence_fails_as_invalid_residue(self, bad_byte_in_seq_file):
        with pytest.raises(ProtfastaException, match='invalid amino acid') as exc:
            protfasta.read_fasta(bad_byte_in_seq_file)
        # the message must itself be printable (no surrogates left in it)
        str(exc.value).encode('utf-8')

    def test_bad_byte_in_sequence_fails_streaming(self, bad_byte_in_seq_file):
        with pytest.raises(ProtfastaException, match='invalid amino acid'):
            list(protfasta.read_fasta_stream(bad_byte_in_seq_file))

    def test_bad_byte_in_sequence_removed(self, bad_byte_in_seq_file):
        assert protfasta.read_fasta(bad_byte_in_seq_file, invalid_sequence_action='remove') == {'h2': 'FGHI'}

    def test_bad_byte_in_sequence_ignored_and_roundtripped(self, bad_byte_in_seq_file, tmp_path):
        result = protfasta.read_fasta(bad_byte_in_seq_file, invalid_sequence_action='ignore')
        assert len(result) == 2
        out = tmp_path / 'out.fasta'
        protfasta.write_fasta(result, out)
        assert b'AC\xe9DE' in out.read_bytes()

    def test_bad_byte_duplicate_messages_are_printable(self, tmp_path):
        f = tmp_path / 'dup.fasta'
        f.write_bytes(b'>h\xe9\nACDE\n>h\xe9\nACDE\n')
        with pytest.raises(ProtfastaException) as exc:
            protfasta.read_fasta(str(f), expect_unique_header=False, duplicate_record_action='fail')
        str(exc.value).encode('utf-8')
        with pytest.raises(ProtfastaException) as exc:
            protfasta.read_fasta(str(f))
        str(exc.value).encode('utf-8')

    def test_utf8_header_preserved(self, tmp_path):
        f = tmp_path / 'utf8.fasta'
        f.write_text('>protéine β\nACDE\n', encoding='utf-8')
        result = protfasta.read_fasta(str(f))
        assert result == {'protéine β': 'ACDE'}
        out = tmp_path / 'out.fasta'
        protfasta.write_fasta(result, out)
        assert out.read_text(encoding='utf-8') == '>protéine β\nACDE\n\n'

    def test_printable_helper(self):
        assert _utilities._printable('plain') == 'plain'
        assert _utilities._printable('é') == 'é'
        assert _utilities._printable('a\udce9b') == 'a\\udce9b'


# ---------------------------------------------------------------------------
# TestParserEdgeCases
# ---------------------------------------------------------------------------
class TestParserEdgeCases:
    """Line-level edge cases in the shared parsing engine."""

    def test_lines_without_newlines(self):
        result = _io._parse_fasta_all(['>h1', 'ACD', 'EFG', '>h2', 'HIK'])
        assert result == [['h1', 'ACDEFG'], ['h2', 'HIK']]

    def test_crlf_lines(self):
        result = _io._parse_fasta_all(['>h1\r\n', 'ACD\r\n', 'EFG\r\n'])
        assert result == [['h1', 'ACDEFG']]

    def test_crlf_file_on_disk(self, tmp_path):
        f = tmp_path / 'crlf.fasta'
        f.write_bytes(b'>h1\r\nACD\r\nEFG\r\n>h2\r\nHIK\r\n')
        assert protfasta.read_fasta(str(f)) == {'h1': 'ACDEFG', 'h2': 'HIK'}
        assert list(protfasta.read_fasta_stream(str(f))) == [('h1', 'ACDEFG'), ('h2', 'HIK')]

    def test_trailing_whitespace_stripped(self):
        result = _io._parse_fasta_all(['>h1   \n', 'ACD  \t\n', 'EFG\n'])
        assert result == [['h1', 'ACDEFG']]

    def test_interior_whitespace_preserved_for_validation(self):
        # interior whitespace is data (an invalid residue), not formatting
        result = _io._parse_fasta_all(['>h1\n', 'AC D\n'])
        assert result == [['h1', 'AC D']]

    def test_whitespace_only_lines_skipped(self):
        result = _io._parse_fasta_all(['>h1\n', '   \n', 'ACD\n', '\t\n', 'EFG\n'])
        assert result == [['h1', 'ACDEFG']]

    def test_leading_junk_before_first_header_ignored(self):
        result = _io._parse_fasta_all(['junk\n', 'more junk\n', '>h1\n', 'ACD\n'])
        assert result == [['h1', 'ACD']]

    def test_header_without_sequence_skipped(self):
        result = _io._parse_fasta_all(['>empty\n', '>h1\n', 'ACD\n', '>trailing_empty\n'])
        assert result == [['h1', 'ACD']]

    def test_lowercase_uppercased(self):
        result = _io._parse_fasta_all(['>h1\n', 'acd\n', 'Efg\n'])
        assert result == [['h1', 'ACDEFG']]

    def test_empty_header_allowed(self):
        result = _io._parse_fasta_all(['>\n', 'ACD\n'])
        assert result == [['', 'ACD']]

    def test_gt_inside_sequence_line_is_data(self):
        result = _io._parse_fasta_all(['>h1\n', 'AC>D\n'])
        assert result == [['h1', 'AC>D']]

    def test_indented_header_is_not_a_header(self):
        result = _io._parse_fasta_all(['>h1\n', 'ACD\n', '  >h2\n', 'EFG\n'])
        assert result == [['h1', 'ACD  >H2EFG']]

    def test_no_trailing_newline_at_eof(self, tmp_path):
        f = tmp_path / 'eof.fasta'
        f.write_bytes(b'>h1\nACD\nEFG')
        assert protfasta.read_fasta(str(f)) == {'h1': 'ACDEFG'}

    def test_duplicate_header_reported_immediately_in_stream(self, tmp_path):
        f = tmp_path / 'dup.fasta'
        f.write_text('>h1\nACD\n>h1\nEFG\n>h2\nHIK\n')
        stream = protfasta.read_fasta_stream(str(f), expect_unique_header=True, silence_warnings=True)
        assert next(stream) == ('h1', 'ACD')
        with pytest.raises(ProtfastaException, match='duplicate header'):
            next(stream)

    def test_iter_fasta_wrapper(self):
        assert list(_io._iter_fasta(SIMPLE_FILE))[0] == (WASL_HEADER, WASL_SEQ)


# ---------------------------------------------------------------------------
# TestConverter
# ---------------------------------------------------------------------------
class TestConverter:
    """The shared converter built by _make_converter."""

    def test_builtin_matches_convert_to_valid(self):
        convert = _utilities._make_converter()
        for seq in ['ACDB', 'A--CD', 'A C D', 'BUXZ*', 'ACDEF']:
            assert convert(seq) == _utilities.convert_to_valid(seq)

    def test_builtin_alignment_matches_convert_to_valid(self):
        convert = _utilities._make_converter(alignment=True)
        for seq in ['ACDB', 'A--CD', 'A C D', 'BUXZ*-', 'ACDEF']:
            assert convert(seq) == _utilities.convert_to_valid(seq, alignment=True)

    def test_custom_single_char(self):
        convert = _utilities._make_converter({'.': 'A', '-': ''})
        assert convert('G.H-I') == 'GAHI'

    def test_custom_multi_char(self):
        convert = _utilities._make_converter({'..': 'A'})
        assert convert('G..H') == 'GAH'
        assert convert('G.H') == 'G.H'

    def test_reusable(self):
        convert = _utilities._make_converter({'.': 'A'})
        assert [convert(s) for s in ['.', '..', 'A']] == ['A', 'AA', 'A']

    def test_convert_invalid_sequences_counts_only_changed(self):
        data = [['h1', 'ACD'], ['h2', 'AC.'], ['h3', '..']]
        _, count = _utilities.convert_invalid_sequences(data, {'.': 'D'})
        assert count == 2
        assert [d[1] for d in data] == ['ACD', 'ACD', 'DD']

    def test_convert_invalid_sequences_identity_mapping_not_counted(self):
        data = [['h1', 'ACD']]
        _, count = _utilities.convert_invalid_sequences(data, {'A': 'A'})
        assert count == 0

    def test_read_fasta_convert_with_multi_char_custom_dict(self, tmp_path):
        f = tmp_path / 'in.fasta'
        f.write_text('>h1\nACXXD\n')
        result = protfasta.read_fasta(str(f), invalid_sequence_action='convert', correction_dictionary={'XX': 'G'})
        assert result == {'h1': 'ACGD'}

    def test_read_fasta_rejects_bad_correction_dictionary_before_reading(self):
        with pytest.raises(ProtfastaException, match='non-empty'):
            protfasta.read_fasta('/definitely/not/a/file.fasta', invalid_sequence_action='convert', correction_dictionary={'': 'A'})


# ---------------------------------------------------------------------------
# TestWriteFastaAtomicity
# ---------------------------------------------------------------------------
class TestWriteFastaAtomicity:
    """write_fasta validates everything before touching the filesystem."""

    def test_empty_sequence_does_not_create_file(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta({'h1': 'ACDE', 'h2': ''}, outfile)
        assert not outfile.exists()

    def test_empty_sequence_does_not_truncate_existing_file(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        protfasta.write_fasta({'keep': 'ACDE'}, outfile)
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta({'h1': 'FGHI', 'h2': ''}, outfile)
        assert protfasta.read_fasta(outfile) == {'keep': 'ACDE'}

    def test_empty_sequence_does_not_partially_append(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        protfasta.write_fasta({'keep': 'ACDE'}, outfile)
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta([['h1', 'FGHI'], ['h2', '']], outfile, append_to_fasta=True)
        assert protfasta.read_fasta(outfile) == {'keep': 'ACDE'}

    def test_bad_list_element_does_not_create_file(self, tmp_path):
        outfile = tmp_path / 'out.fasta'
        with pytest.raises(ProtfastaException):
            protfasta.write_fasta([['h1', 'ACDE'], ['h2']], outfile)
        assert not outfile.exists()

    def test_unwritable_path_is_protfasta_exception(self, tmp_path):
        with pytest.raises(ProtfastaException, match='Unable to open file for writing'):
            protfasta.write_fasta({'h1': 'ACDE'}, tmp_path / 'no_such_dir' / 'out.fasta')

    def test_read_fasta_unwritable_output_is_protfasta_exception(self, tmp_path):
        with pytest.raises(ProtfastaException, match='Unable to open file for writing'):
            protfasta.read_fasta(SIMPLE_FILE, output_filename=tmp_path / 'no_such_dir' / 'out.fasta')


# ---------------------------------------------------------------------------
# TestStreamRobustness
# ---------------------------------------------------------------------------
class TestStreamRobustness:
    """File-handling edge cases specific to read_fasta_stream."""

    def test_missing_file_raises_at_call_time(self):
        # no need to iterate - a missing input is an argument error
        with pytest.raises(ProtfastaException, match='Unable to find file'):
            protfasta.read_fasta_stream('does_not_exist_12345.fasta')

    def test_missing_input_does_not_create_output(self, tmp_path):
        out = tmp_path / 'out.fasta'
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream('does_not_exist_12345.fasta', output_filename=out))
        assert not out.exists()

    def test_directory_input_does_not_create_output(self, tmp_path):
        out = tmp_path / 'out.fasta'
        with pytest.raises(ProtfastaException):
            list(protfasta.read_fasta_stream(str(tmp_path), output_filename=out))
        assert not out.exists()

    def test_unwritable_output_is_protfasta_exception(self, tmp_path):
        with pytest.raises(ProtfastaException, match='Unable to open file for writing'):
            list(protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=tmp_path / 'no_such_dir' / 'out.fasta'))

    def test_symlinked_output_same_as_input_raises(self, tmp_path):
        link = tmp_path / 'link.fasta'
        try:
            os.symlink(SIMPLE_FILE, link)
        except (OSError, NotImplementedError):
            pytest.skip('symlinks not supported here')
        with pytest.raises(ProtfastaException, match='must differ'):
            protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=link)

    def test_tee_output_identical_to_write_fasta(self, tmp_path):
        a = tmp_path / 'stream.fasta'
        b = tmp_path / 'write.fasta'
        records = list(protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=a, return_list=True))
        protfasta.write_fasta(records, b)
        assert a.read_bytes() == b.read_bytes()

    def test_input_handle_closed_after_exhaustion(self, tmp_path):
        import gc
        stream = protfasta.read_fasta_stream(SIMPLE_FILE)
        list(stream)
        gc.collect()
        # generator is finished; closing again is a no-op and must not raise
        stream.close()

    def test_partial_consumption_then_close_completes_output_so_far(self, tmp_path):
        out = tmp_path / 'partial.fasta'
        stream = protfasta.read_fasta_stream(SIMPLE_FILE, output_filename=out)
        first = next(stream)
        stream.close()
        assert protfasta.read_fasta(out, return_list=True) == [list(first)]


# ---------------------------------------------------------------------------
# TestReadFastaStreamParity
# ---------------------------------------------------------------------------
class TestReadFastaStreamParity:
    """read_fasta and read_fasta_stream must agree wherever both succeed."""

    FILES = [
        SIMPLE_FILE, DUPLICATE_RECORD_FILE, DUPLICATE_SEQ_FILE, BADCHAR_FILE,
        NONSTANDARD_FILE, FIXABLE_INVALID_FILE, UNFIXABLE_INVALID_FILE,
        ALIGNED_VALID_FILE, ALIGNED_CONVERTABLE_FILE, ALIGNED_UNCONVERTABLE_FILE,
    ]
    ACTIONS = ['ignore', 'fail', 'remove', 'convert', 'convert-ignore', 'convert-remove']

    @pytest.mark.parametrize('alignment', [False, True])
    @pytest.mark.parametrize('action', ACTIONS)
    def test_parity_over_all_test_files(self, action, alignment):
        for fn in self.FILES:
            kwargs = dict(
                expect_unique_header=False,
                duplicate_record_action='remove',
                duplicate_sequence_action='remove',
                invalid_sequence_action=action,
                alignment=alignment,
                return_list=True,
            )
            try:
                ref = protfasta.read_fasta(fn, **kwargs)
            except ProtfastaException:
                with pytest.raises(ProtfastaException):
                    list(protfasta.read_fasta_stream(fn, silence_warnings=True, **kwargs))
                continue
            streamed = list(protfasta.read_fasta_stream(fn, silence_warnings=True, **kwargs))
            assert streamed == ref, fn
