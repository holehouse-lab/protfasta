"""
Internal utility functions for sequence validation, conversion, and
duplicate handling used by the protfasta processing pipeline.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""

from __future__ import annotations

import hashlib
from typing import Callable, Mapping, Optional, Union

from .protfasta_exceptions import ProtfastaException
from ._configs import (
    STANDARD_CONVERSION,
    STANDARD_CONVERSION_WITH_GAP,
    _TRANSLATE_STANDARD,
    _TRANSLATE_WITH_GAP,
    _VALID_BYTES_STANDARD,
    _VALID_BYTES_WITH_GAP,
    _DELETE_VALID_STANDARD,
    _DELETE_VALID_WITH_GAP,
)


####################################################################################################
#
#
def _seq_hash(seq: str) -> bytes:
    """Return a 16-byte blake2b digest of *seq* for cheap duplicate lookup.

    Used by duplicate-detection utilities so the lookup structure stores
    128-bit digests instead of whole (potentially very long) sequences.
    Collision probability for 10**8 sequences is ~10**-22.

    Note that sequences are encoded as UTF-8 rather than ASCII.  Duplicate
    detection runs *before* invalid-residue handling, so it must be able to
    digest whatever the file contained -- including non-ASCII junk that a
    later ``invalid_sequence_action`` will flag or strip.  The
    ``surrogatepass`` handler additionally lets lone surrogates (which is
    how undecodable bytes arrive from a file read with ``surrogateescape``)
    through unharmed.

    Parameters
    ----------
    seq : str
        Amino acid sequence to digest.

    Returns
    -------
    bytes
        A 16-byte digest.
    """
    return hashlib.blake2b(seq.encode('utf-8', 'surrogatepass'), digest_size=16).digest()


####################################################################################################
#
#
def _record_hash(header: str, seq: str) -> bytes:
    """Return a 16-byte blake2b digest identifying a (header, sequence) record.

    Duplicate *record* detection needs a key that combines the header and
    the sequence, so that a header seen with several different sequences
    is still caught if any one of those sequences is later repeated.  The
    header and sequence are fed through a single hash with a NUL separator
    (a header can never contain one, since it is a single line of text
    with trailing whitespace stripped) so the cost is one pass over the
    record, the same as :func:`_seq_hash`.

    Parameters
    ----------
    header : str
        The record header.

    seq : str
        The sequence.

    Returns
    -------
    bytes
        A 16-byte digest unique to the (header, sequence) pair.
    """
    h = hashlib.blake2b(digest_size=16)
    h.update(header.encode('utf-8', 'surrogatepass'))
    h.update(b'\0')
    h.update(seq.encode('utf-8', 'surrogatepass'))
    return h.digest()


####################################################################################################
#
#
def _printable(s: str) -> str:
    """Return *s* in a form that is safe to embed in an error message.

    Files are read with the ``surrogateescape`` error handler, so bytes
    that were not valid UTF-8 survive as lone surrogate code points.  Those
    cannot be printed to a UTF-8 terminal, so a message that embedded one
    verbatim would itself raise a ``UnicodeEncodeError`` the moment a user
    tried to print it.  This converts any such code points into
    ``\\udcXX`` escapes and leaves everything else untouched.

    Parameters
    ----------
    s : str
        Text destined for an exception message.

    Returns
    -------
    str
        Printable text.
    """
    return s.encode('utf-8', 'backslashreplace').decode('utf-8')


####################################################################################################
#
#
def _validate_correction_dictionary(correction_dictionary: dict[str, str]) -> None:
    """Check that a correction dictionary maps non-empty strings to strings.

    Anything else either crashes ``str.maketrans`` with an opaque error or,
    worse, silently corrupts sequences: an empty-string key makes
    ``str.replace`` insert the replacement between every pair of residues.

    Parameters
    ----------
    correction_dictionary : dict
        The user-supplied mapping.

    Raises
    ------
    ProtfastaException
        If any key is not a non-empty string or any value is not a string.
    """
    for key, value in correction_dictionary.items():
        if not isinstance(key, str) or len(key) == 0:
            raise ProtfastaException("keys of 'correction_dictionary' must be non-empty strings (got %r)" % (key,))
        if not isinstance(value, str):
            raise ProtfastaException("values of 'correction_dictionary' must be strings (got %r for key %r)" % (value, key))


####################################################################################################
#
#
def build_custom_dictionary(additional_dictionary: dict[str, str]) -> dict[str, str]:
    """Build a correction dictionary by merging defaults with custom entries.

    Starts from the built-in ``STANDARD_CONVERSION`` table and overlays
    *additional_dictionary* on top, so caller-supplied mappings take
    precedence over the defaults.

    Parameters
    ----------
    additional_dictionary : dict[str, str]
        Mapping of non-standard characters to their desired
        replacements.  Keys already present in the default table
        will be overwritten.

    Returns
    -------
    dict[str, str]
        The merged correction dictionary.
    """
    # the additional dictionary comes second so it deliberately overwrites
    # the standard conversion, i.e. a passed dictionary takes precedence
    return {**STANDARD_CONVERSION, **additional_dictionary}


####################################################################################################
#
#
def _make_converter(
    correction_dictionary: Optional[dict[str, str]] = None,
    alignment: bool = False,
) -> Callable[[str], str]:
    """Build a function that converts non-standard residues in a sequence.

    The returned callable applies *correction_dictionary* (or the built-in
    table when that is ``None``/empty) to a sequence.  Building it once and
    reusing it across every record avoids rebuilding a ``str.translate``
    table per sequence, and the converter itself is guarded: it first
    checks whether any convertible character is actually present (a handful
    of C-level ``in`` scans) and returns the input object unchanged if not.
    On typical data, where the vast majority of sequences are already
    clean, that skips the translate pass and its allocation entirely.

    Parameters
    ----------
    correction_dictionary : dict[str, str] or None, optional
        Mapping of characters (or multi-character strings) to their
        replacements.  ``None`` or an empty dictionary selects the
        built-in table.

    alignment : bool, optional
        When ``True`` the built-in table preserves dashes (``'-'``) as
        gap characters.  Ignored when a custom dictionary is supplied.
        Default ``False``.

    Returns
    -------
    callable
        A function ``(str) -> str``.

    Raises
    ------
    ProtfastaException
        If *correction_dictionary* contains an empty key or a non-string
        key or value.
    """

    table: Optional[Mapping[int, Union[str, int, None]]]
    if not correction_dictionary:
        if alignment:
            table = _TRANSLATE_WITH_GAP
            keys = tuple(STANDARD_CONVERSION_WITH_GAP)
        else:
            table = _TRANSLATE_STANDARD
            keys = tuple(STANDARD_CONVERSION)
    else:
        _validate_correction_dictionary(correction_dictionary)
        keys = tuple(correction_dictionary)

        # str.translate can only map single characters (its table is keyed
        # by code point); a dictionary with any multi-character key falls
        # back to sequential str.replace calls, applied in dictionary order
        # as it always has been.
        if all(len(k) == 1 for k in keys):
            table = {ord(k): v for k, v in correction_dictionary.items()}
        else:
            table = None

    if table is not None:
        tbl = table

        def convert(seq: str) -> str:
            for k in keys:
                if k in seq:
                    return seq.translate(tbl)
            return seq

    else:
        replacements = tuple(correction_dictionary.items()) if correction_dictionary else ()

        def convert(seq: str) -> str:
            for k in keys:
                if k in seq:
                    for old, new in replacements:
                        seq = seq.replace(old, new)
                    return seq
            return seq

    return convert


# Converters for the two built-in tables, built once at import time so
# convert_to_valid() pays nothing for the common no-custom-dictionary case.
_CONVERT_STANDARD = _make_converter(None, alignment=False)
_CONVERT_WITH_GAP = _make_converter(None, alignment=True)


####################################################################################################
#
#
def convert_to_valid(
    seq: str,
    correction_dictionary: Optional[dict[str, str]] = None,
    alignment: bool = False,
) -> str:
    """Convert non-standard amino-acid characters to standard ones.

    Default conversions (when no *correction_dictionary* is supplied):

    * ``B`` -> ``N``
    * ``U`` -> ``C``
    * ``X`` -> ``G``
    * ``Z`` -> ``Q``
    * ``' '`` -> ``''`` (space removed)
    * ``*`` -> ``''``
    * ``-`` -> ``''`` (only when *alignment* is ``False``)

    When converting many sequences with the same custom dictionary prefer
    :func:`convert_invalid_sequences` (or :func:`_make_converter`
    directly), which builds the conversion table once rather than per
    call.

    Parameters
    ----------
    seq : str
        Amino acid sequence to convert.

    correction_dictionary : dict[str, str] or None, optional
        Custom mapping of characters to replacements.  When provided
        this is used instead of the built-in table.

    alignment : bool, optional
        When ``True``, dashes (``'-'``) are kept as valid gap
        characters.  Default ``False``.

    Returns
    -------
    str
        The sequence with non-standard residues replaced.  If nothing
        needed converting the input object itself is returned.

    Raises
    ------
    ProtfastaException
        If *correction_dictionary* contains an empty key or a non-string
        key or value.
    """
    if not correction_dictionary:
        return _CONVERT_WITH_GAP(seq) if alignment else _CONVERT_STANDARD(seq)
    return _make_converter(correction_dictionary, alignment)(seq)


####################################################################################################
#
#
def check_sequence_is_valid(
    seq: str,
    alignment: bool = False,
) -> tuple[bool, Union[str, int]]:
    """Check whether every character in *seq* is a standard amino acid.

    The check deletes every valid residue from the sequence in a single
    C-level pass; whatever is left over is, by construction, the list of
    invalid characters in order of appearance.  Pure-ASCII sequences (the
    overwhelmingly common case) go through ``bytes.translate``, which is
    several times faster again than the ``str`` equivalent.

    Parameters
    ----------
    seq : str
        Amino acid sequence to validate.

    alignment : bool, optional
        When ``True``, dashes (``'-'``) are accepted as valid gap
        characters.  Default ``False``.

    Returns
    -------
    tuple[bool, str | int]
        A two-element tuple:

        * ``(True, 0)`` if the sequence is entirely valid.
        * ``(False, <char>)`` where ``<char>`` is the first invalid
          character encountered.
    """
    if seq.isascii():
        bad = seq.encode('ascii').translate(None, _VALID_BYTES_WITH_GAP if alignment else _VALID_BYTES_STANDARD)
        if bad:
            return (False, chr(bad[0]))
        return (True, 0)

    bad_str = seq.translate(_DELETE_VALID_WITH_GAP if alignment else _DELETE_VALID_STANDARD)
    if bad_str:
        return (False, bad_str[0])
    return (True, 0)


####################################################################################################
#
#
def convert_invalid_sequences(
    dataset: list[list[str]],
    correction_dictionary: Optional[dict[str, str]] = None,
    alignment: bool = False,
) -> tuple[list[list[str]], int]:
    """Convert invalid residues in every sequence in *dataset*.

    Each sequence is passed through a converter built once by
    :func:`_make_converter`.  The dataset is modified in place and also
    returned.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    correction_dictionary : dict[str, str] or None, optional
        Custom character-replacement mapping.  ``None`` uses the
        built-in default table.

    alignment : bool, optional
        When ``True``, dashes are preserved.  Default ``False``.

    Returns
    -------
    tuple[list[list[str]], int]
        A tuple of the (mutated) dataset and the number of sequences
        that were altered.
    """
    convert = _make_converter(correction_dictionary, alignment)

    count = 0
    for entry in dataset:
        s = entry[1]
        new = convert(s)
        # the converter hands back the very same object when it had
        # nothing to do, so the identity test short-circuits the common
        # case and the equality test covers a no-op mapping (e.g. A->A)
        if new is not s and new != s:
            entry[1] = new
            count += 1

    return (dataset, count)


####################################################################################################
#
#
def remove_invalid_sequences(
    dataset: list[list[str]],
    alignment: bool = False,
) -> list[list[str]]:
    """Return only entries whose sequences are fully valid.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    alignment : bool, optional
        When ``True``, dashes are treated as valid.  Default ``False``.

    Returns
    -------
    list[list[str]]
        Filtered list containing only entries with valid sequences.
    """
    return [element for element in dataset if check_sequence_is_valid(element[1], alignment)[0]]


####################################################################################################
#
#
def fail_on_invalid_sequences(
    dataset: list[list[str]],
    alignment: bool = False,
) -> None:
    """Raise if any sequence in *dataset* contains invalid residues.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    alignment : bool, optional
        When ``True``, dashes are treated as valid.  Default ``False``.

    Raises
    ------
    ProtfastaException
        On the first sequence that contains an invalid character.
    """
    for entry in dataset:
        (status, info) = check_sequence_is_valid(entry[1], alignment)
        if status is not True:
            raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (_printable(str(info)), _printable(entry[0]), _printable(entry[1])))


####################################################################################################
#
#
def convert_list_to_dictionary(
    raw_list: list[list[str]],
    verbose: bool = False,
) -> dict[str, str]:
    """Convert a list of ``[header, sequence]`` pairs into a dictionary.

    When *verbose* is ``True``, warnings are printed for overwritten
    duplicate headers and a summary line is emitted.

    Parameters
    ----------
    raw_list : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    verbose : bool, optional
        If ``True``, print duplicate-header warnings to stdout.
        Default ``False``.

    Returns
    -------
    dict[str, str]
        Mapping of header to sequence.  If duplicate headers exist,
        the last occurrence wins.
    """
    if not verbose:
        # dict() consumes an iterable of pairs in C; last occurrence wins,
        # exactly as the explicit loop below
        return dict(raw_list)

    return_dict: dict[str, str] = {}
    warning_count = 0
    for entry in raw_list:
        if entry[0] in return_dict:
            warning_count = warning_count + 1
            print('[WARNING]: Overwriting entry [count = %i]' % (warning_count))
        return_dict[entry[0]] = entry[1]
    if warning_count > 0:
        print('[INFO] If you want to avoid overwriting duplicate headers set return_list=True')
    else:
        print('[INFO]: All processed sequences uniquely added to the returning dictionary')

    return return_dict


####################################################################################################
#
#
def fail_on_duplicates(dataset: list[list[str]]) -> None:
    """Raise if any exact duplicate record exists in *dataset*.

    A duplicate record is defined as two entries with the same header
    **and** the same sequence.  Entries that share a header but have
    different sequences are *not* considered duplicates -- and, crucially,
    a header that appears with several different sequences must still be
    caught if any one of those sequences is later repeated.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        On the first duplicate record found.
    """
    # Each record is reduced to a single 16-byte digest of its header and
    # sequence, so the bookkeeping is one small bytes object per record
    # regardless of how long the sequences are.
    seen: set[bytes] = set()
    for entry in dataset:
        key = _record_hash(entry[0], entry[1])
        if key in seen:
            raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (_printable(entry[0]), _printable(entry[1])))
        seen.add(key)


####################################################################################################
#
#
def remove_duplicates(dataset: list[list[str]]) -> list[list[str]]:
    """Remove exact duplicate records, keeping the first occurrence.

    A duplicate record is defined as two entries with the same header
    **and** the same sequence.  Entries that share a header but have
    different sequences are *not* considered duplicates.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Returns
    -------
    list[list[str]]
        De-duplicated list, preserving original order.
    """
    # One 16-byte record digest per entry (see fail_on_duplicates) rather
    # than a per-header set of sequence digests -- roughly a quarter of the
    # memory on files with hundreds of millions of records.
    seen: set[bytes] = set()
    updated: list[list[str]] = []

    for entry in dataset:
        key = _record_hash(entry[0], entry[1])
        if key not in seen:
            seen.add(key)
            updated.append(entry)

    return updated


####################################################################################################
#
#
def fail_on_duplicate_sequences(dataset: list[list[str]]) -> None:
    """Raise if any two entries share the same sequence.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        On the first pair of entries that share a sequence.
    """
    # Key by a 16-byte digest rather than the full sequence for memory
    # efficiency on large files.
    seq_to_header: dict[bytes, str] = {}
    for entry in dataset:
        digest = _seq_hash(entry[1])
        if digest in seq_to_header:
            raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (_printable(seq_to_header[digest]), _printable(entry[0])))
        seq_to_header[digest] = entry[0]


####################################################################################################
#
#
def remove_duplicate_sequences(dataset: list[list[str]]) -> list[list[str]]:
    """Remove entries with duplicate sequences, keeping the first occurrence.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Returns
    -------
    list[list[str]]
        Filtered list with unique sequences, preserving original order.
    """
    # Track seen sequences by their 16-byte digests to keep peak memory
    # low for files with long sequences.
    lookup: set[bytes] = set()
    updated: list[list[str]] = []

    for entry in dataset:
        digest = _seq_hash(entry[1])
        if digest in lookup:
            continue
        lookup.add(digest)
        updated.append(entry)
    return updated
