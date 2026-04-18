"""
Internal utility functions for sequence validation, conversion, and
duplicate handling used by the protfasta processing pipeline.
"""

from __future__ import annotations

from typing import Optional

import hashlib

from .protfasta_exceptions import ProtfastaException
from ._configs import (
    STANDARD_CONVERSION,
    STANDARD_CONVERSION_WITH_GAP,
    STANDARD_AAS,
    STANDARD_AAS_WITH_GAP,
    _STANDARD_AAS_SET,
    _STANDARD_AAS_WITH_GAP_SET,
    _TRANSLATE_STANDARD,
    _TRANSLATE_WITH_GAP,
)


def _seq_hash(seq: str) -> bytes:
    """Return a 16-byte blake2b digest of *seq* for cheap duplicate lookup.

    Used by duplicate-detection utilities so the lookup structure stores
    128-bit digests instead of whole (potentially very long) sequences.
    Collision probability for 10**8 sequences is ~10**-22.
    """
    return hashlib.blake2b(seq.encode('ascii', 'surrogatepass'), digest_size=16).digest()






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

    final_dict = {}
    for i in STANDARD_CONVERSION:        
        final_dict[i] = STANDARD_CONVERSION[i]

    # this deliberatly overwirtes the standard conversion
    # so a passed additional dictionary takes precedence! 
    for i in additional_dictionary:
        final_dict[i] = additional_dictionary[i]
        
    return final_dict



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
        The sequence with non-standard residues replaced.
    """

    # Fast path: use the pre-built str.translate tables for the built-in
    # conversion dictionaries.  str.translate is a single C-level pass
    # and is typically 5-20x faster than chained str.replace calls.
    if not correction_dictionary:
        if alignment:
            return seq.translate(_TRANSLATE_WITH_GAP)
        return seq.translate(_TRANSLATE_STANDARD)

    # Custom dictionary path: build a translate table on the fly if every
    # key is a single character (the common case), else fall back to the
    # original str.replace loop which handles multi-character keys.
    if all(len(k) == 1 for k in correction_dictionary):
        return seq.translate(str.maketrans(correction_dictionary))

    for i in correction_dictionary:
        seq = seq.replace(i, correction_dictionary[i])
    return seq



####################################################################################################
#
#    
def check_sequence_is_valid(
    seq: str,
    alignment: bool = False,
) -> tuple[bool, str | int]:
    """Check whether every character in *seq* is a standard amino acid.

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

    valid = _STANDARD_AAS_WITH_GAP_SET if alignment else _STANDARD_AAS_SET

    # Compute the set of distinct invalid characters via set difference
    # (a single C-level pass to build set(seq), then frozenset-diff).
    invalid = set(seq).difference(valid)
    if invalid:
        # Return a deterministic first-invalid character: scan seq once
        # so that the reported offender matches the original ordering.
        for ch in seq:
            if ch in invalid:
                return (False, ch)

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

    Each sequence is passed through :func:`convert_to_valid`.  The
    dataset is modified in place and also returned.

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

    count = 0
    for idx in range(0,len(dataset)):
        s = dataset[idx][1]
        dataset[idx][1] = convert_to_valid(s, correction_dictionary, alignment)
        
        if s != dataset[idx][1]:
            count = count + 1
        
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
            raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (info, entry[0], entry[1]))



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
    if verbose:
        return_dict={}
        warning_count=0
        for entry in raw_list:
            if entry[0] in return_dict:
                warning_count=warning_count+1
                print('[WARNING]: Overwriting entry [count = %i]'%(warning_count))
            return_dict[entry[0]] = entry[1]
        if warning_count > 0:
            print('[INFO] If you want to avoid overwriting duplicate headers set return_list=True')
        else:
            print('[INFO]: All processed sequences uniquely added to the returning dictionary')
            
    else:
        return_dict={}
        for entry in raw_list:
            return_dict[entry[0]] = entry[1]

    return return_dict



####################################################################################################
#
#    
def fail_on_duplicates(dataset: list[list[str]]) -> None:
    """Raise if any exact duplicate record exists in *dataset*.

    A duplicate record is defined as two entries with the same header
    **and** the same sequence.

    Parameters
    ----------
    dataset : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        On the first duplicate record found.
    """
    # Store blake2b digests of sequences instead of full sequences to
    # keep peak memory low for files with very long sequences.
    lookup: dict[str, bytes] = {}
    for entry in dataset:
        digest = _seq_hash(entry[1])
        if entry[0] not in lookup:
            lookup[entry[0]] = digest
        else:
            if lookup[entry[0]] == digest:
                raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (entry[0], entry[1]))



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
    # Store sets of sequence hashes per header (instead of raw sequences)
    # to keep peak memory low.  O(1) membership test per header.
    lookup: dict[str, set[bytes]] = {}
    updated: list[list[str]] = []

    for entry in dataset:
        header = entry[0]
        digest = _seq_hash(entry[1])

        seen = lookup.get(header)
        if seen is None:
            lookup[header] = {digest}
            updated.append(entry)
        elif digest not in seen:
            seen.add(digest)
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
            raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (seq_to_header[digest], entry[0]))
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
            

    
                

                
