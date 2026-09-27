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

    Used by the streaming duplicate checks in :func:`io._stream_fasta`, which
    must remember every sequence they have seen without keeping the
    sequences themselves (that would defeat the point of streaming).  The
    load-everything path does not need it: there every sequence is already
    in memory, so :func:`fail_on_duplicate_sequences` and
    :func:`remove_duplicate_sequences` compare the strings directly.
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

    Used by the streaming duplicate-record check in
    :func:`io._stream_fasta` (the load-everything path compares records
    directly; see :func:`fail_on_duplicates`).  The key has to combine the
    header and the sequence, so that a header seen with several different
    sequences is still caught if any one of those sequences is later
    repeated, and the two are fed through a single hash so the cost is one
    pass over the record, the same as :func:`_seq_hash`.

    The encoded header is preceded by its length.  Without that the split
    between header and sequence would be ambiguous whatever separator was
    used, because a header can contain any character a file can hold -
    with a NUL separator, for instance, the records ``('A\\0', 'B')`` and
    ``('A', '\\0B')`` hashed identically and one of them was discarded as a
    duplicate of the other.

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
    header_bytes = header.encode('utf-8', 'surrogatepass')
    h = hashlib.blake2b(len(header_bytes).to_bytes(8, 'little'), digest_size=16)
    h.update(header_bytes)
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


# Number of records handled per batch by the dataset-level validation and
# conversion helpers below. Each batch is joined into one string and
# examined with a handful of C-level passes, so the Python-level overhead
# of a per-sequence call is paid once per batch rather than once per record
# whenever the batch turns out to be clean (the overwhelmingly common
# case). The value is not critical - anything from a few hundred to a few
# thousand measures the same - and the transient joined string is only a
# couple of hundred kilobytes at 512 records of typical length.
_BATCH = 512


####################################################################################################
#
#
def _contains_any(text: str, keys: tuple[str, ...]) -> bool:
    """Return ``True`` if any of *keys* occurs in *text*.

    A loop of ``in`` tests, each a single C-level search, which is far
    cheaper than any Python-level scan over *text*.

    Parameters
    ----------
    text : str
        The string to search.

    keys : tuple[str, ...]
        The substrings to look for.

    Returns
    -------
    bool
        Whether at least one key is present.
    """
    for k in keys:
        if k in text:
            return True
    return False


####################################################################################################
#
#
def _batch_is_valid(batch: list[list[str]], valid_bytes: bytes) -> bool:
    """Return ``True`` if every sequence in *batch* is entirely valid.

    The sequences are concatenated and checked in one pass, exactly as
    :func:`check_sequence_is_valid` checks a single sequence.  No separator
    is needed: deleting every valid residue from the concatenation leaves
    behind precisely the invalid characters of every member, so nothing can
    hide across a boundary between two sequences.  A ``False`` result only
    says that at least one member is invalid; the caller falls back to
    per-sequence checks to find out which.

    Parameters
    ----------
    batch : list[list[str]]
        A slice of parsed FASTA data -- ``[header, sequence]`` pairs.

    valid_bytes : bytes
        The ASCII alphabet of valid residues (one of the
        ``_VALID_BYTES_*`` tables).

    Returns
    -------
    bool
        Whether every sequence in the batch is valid.
    """
    joined = ''.join([entry[1] for entry in batch])
    return joined.isascii() and not joined.encode('ascii').translate(None, valid_bytes)


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
def _conversion_keys(
    correction_dictionary: Optional[dict[str, str]] = None,
    alignment: bool = False,
) -> tuple[str, ...]:
    """Return the characters (or strings) a converter would replace.

    This is the set of keys of the table that :func:`_make_converter` would
    apply for the same arguments: the caller's *correction_dictionary* when
    one is given, otherwise the built-in table (with or without the gap
    character, depending on *alignment*).

    Parameters
    ----------
    correction_dictionary : dict[str, str] or None, optional
        Custom mapping.  ``None`` or empty selects the built-in table.

    alignment : bool, optional
        Whether the built-in table should preserve dashes.  Ignored when a
        custom dictionary is supplied.  Default ``False``.

    Returns
    -------
    tuple[str, ...]
        The keys, in table order.
    """
    if not correction_dictionary:
        return tuple(STANDARD_CONVERSION_WITH_GAP if alignment else STANDARD_CONVERSION)
    return tuple(correction_dictionary)


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
        table = _TRANSLATE_WITH_GAP if alignment else _TRANSLATE_STANDARD
    else:
        _validate_correction_dictionary(correction_dictionary)

        # str.translate can only map single characters (its table is keyed
        # by code point); a dictionary with any multi-character key falls
        # back to sequential str.replace calls, applied in dictionary order
        # as it always has been.
        if all(len(k) == 1 for k in correction_dictionary):
            table = {ord(k): v for k, v in correction_dictionary.items()}
        else:
            table = None

    keys = _conversion_keys(correction_dictionary, alignment)

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

    Sequences are processed in batches: the members of a batch are joined
    and scanned once for each convertible character, and a batch in which
    none occurs (the common case on real data) is skipped without touching
    any of its sequences individually.

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
    keys = _conversion_keys(correction_dictionary, alignment)

    count = 0
    for start in range(0, len(dataset), _BATCH):
        batch = dataset[start:start + _BATCH]

        # A multi-character key can straddle the boundary between two
        # adjacent sequences in the joined string, which only ever makes
        # this test err on the side of looking at the batch member by
        # member - never on the side of skipping a sequence that needs
        # converting.
        if not _contains_any(''.join([entry[1] for entry in batch]), keys):
            continue

        for entry in batch:
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

    Sequences are checked in batches (see :func:`_batch_is_valid`); only a
    batch that contains at least one invalid sequence is examined member by
    member.

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
    valid_bytes = _VALID_BYTES_WITH_GAP if alignment else _VALID_BYTES_STANDARD

    updated: list[list[str]] = []
    for start in range(0, len(dataset), _BATCH):
        batch = dataset[start:start + _BATCH]
        if _batch_is_valid(batch, valid_bytes):
            updated.extend(batch)
        else:
            updated.extend([entry for entry in batch if check_sequence_is_valid(entry[1], alignment)[0]])

    return updated


####################################################################################################
#
#
def fail_on_invalid_sequences(
    dataset: list[list[str]],
    alignment: bool = False,
) -> None:
    """Raise if any sequence in *dataset* contains invalid residues.

    Sequences are checked in batches (see :func:`_batch_is_valid`); only a
    batch that contains at least one invalid sequence is examined member by
    member, so the sequence reported is always the first invalid one in
    dataset order.

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
    valid_bytes = _VALID_BYTES_WITH_GAP if alignment else _VALID_BYTES_STANDARD

    for start in range(0, len(dataset), _BATCH):
        batch = dataset[start:start + _BATCH]
        if _batch_is_valid(batch, valid_bytes):
            continue
        for entry in batch:
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
def _seen_record(seen: dict[str, Union[str, set[str]]], header: str, seq: str) -> bool:
    """Note the record ``(header, seq)`` in *seen*, reporting whether it was already there.

    This is the bookkeeping behind :func:`fail_on_duplicates` and
    :func:`remove_duplicates`.  Every sequence is already in memory on this
    path, so records are compared directly - exactly, with no digests - and
    the bookkeeping only ever holds references to strings the dataset
    already owns.

    *seen* is keyed by header.  Its value is the first sequence seen under
    that header, and only if the header turns up again with a *different*
    sequence is the value promoted to a set of every sequence seen under it.
    On typical data, where a header recurs rarely if at all, that means one
    dictionary entry per header, only the header is ever hashed, and a
    repeated header is checked with a single string comparison.

    Parameters
    ----------
    seen : dict[str, str or set[str]]
        Bookkeeping built up by previous calls; updated in place.

    header : str
        The record header.

    seq : str
        The record sequence.

    Returns
    -------
    bool
        ``True`` if this exact record has been seen before, otherwise
        ``False`` (in which case it has now been noted).
    """
    prior = seen.get(header)
    if prior is None:
        seen[header] = seq
        return False
    if isinstance(prior, str):
        if prior == seq:
            return True
        seen[header] = {prior, seq}
        return False
    if seq in prior:
        return True
    prior.add(seq)
    return False


####################################################################################################
#
#
def fail_on_duplicates(dataset: list[list[str]]) -> None:
    """Raise if any exact duplicate record exists in *dataset*.

    A duplicate record is defined as two entries with the same header
    **and** the same sequence.  Entries that share a header but have
    different sequences are *not* considered duplicates -- and, crucially,
    a header that appears with several different sequences must still be
    caught if any one of those sequences is later repeated.  Records are
    compared exactly (see :func:`_seen_record`).

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        On the first duplicate record found.
    """
    seen: dict[str, Union[str, set[str]]] = {}
    for entry in dataset:
        if _seen_record(seen, entry[0], entry[1]):
            raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (_printable(entry[0]), _printable(entry[1])))


####################################################################################################
#
#
def remove_duplicates(dataset: list[list[str]]) -> list[list[str]]:
    """Remove exact duplicate records, keeping the first occurrence.

    A duplicate record is defined as two entries with the same header
    **and** the same sequence.  Entries that share a header but have
    different sequences are *not* considered duplicates.  Records are
    compared exactly (see :func:`_seen_record`).

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Returns
    -------
    list[list[str]]
        De-duplicated list, preserving original order.
    """
    seen: dict[str, Union[str, set[str]]] = {}
    return [entry for entry in dataset if not _seen_record(seen, entry[0], entry[1])]


####################################################################################################
#
#
def fail_on_duplicate_sequences(dataset: list[list[str]]) -> None:
    """Raise if any two entries share the same sequence.

    Sequences are compared exactly, by keying a dictionary on the sequence
    strings the dataset already holds (so no sequence is copied, and each
    one is hashed once, by Python's own string hash).

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        On the first pair of entries that share a sequence.
    """
    # sequence -> the first header it was seen under, so the error can name
    # both records
    first_header: dict[str, str] = {}
    for entry in dataset:
        seq = entry[1]
        if seq in first_header:
            raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (_printable(first_header[seq]), _printable(entry[0])))
        first_header[seq] = entry[0]


####################################################################################################
#
#
def remove_duplicate_sequences(dataset: list[list[str]]) -> list[list[str]]:
    """Remove entries with duplicate sequences, keeping the first occurrence.

    Sequences are compared exactly, through a set of the sequence strings
    the dataset already holds (so no sequence is copied).

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Returns
    -------
    list[list[str]]
        Filtered list with unique sequences, preserving original order.
    """
    seen: set[str] = set()
    add = seen.add
    updated: list[list[str]] = []

    for entry in dataset:
        seq = entry[1]
        if seq in seen:
            continue
        add(seq)
        updated.append(entry)
    return updated
