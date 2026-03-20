"""
Internal utility functions for sequence validation, conversion, and
duplicate handling used by the protfasta processing pipeline.
"""

from __future__ import annotations

from typing import Optional

from .protfasta_exceptions import ProtfastaException
from ._configs import STANDARD_CONVERSION, STANDARD_CONVERSION_WITH_GAP, STANDARD_AAS, STANDARD_AAS_WITH_GAP






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

    if correction_dictionary:
        converter = correction_dictionary
    else:
        if alignment:
            converter = STANDARD_CONVERSION_WITH_GAP
        else:
            converter = STANDARD_CONVERSION

    # cycle over each key in the converter dictionary and replace all values in
    # the sequence with the corresponding converted one 
    for i in converter:
        seq = seq.replace('%s'%(i), converter[i])

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

    if alignment == True:
        valid_AA_list = STANDARD_AAS_WITH_GAP
    else:
        valid_AA_list = STANDARD_AAS
    
    s = list(set(seq))
    for i in s:
        if i not in valid_AA_list:
            return (False, i)

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
    lookup: dict[str, str] = {}
    for entry in dataset:
        if entry[0] not in lookup:
            lookup[entry[0]] = entry[1]
        else:
            if lookup[entry[0]] == entry[1]:
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
    lookup: dict[str, list[str]] = {}
    updated: list[list[str]] = []

    # for each entry in the dataset
    for entry in dataset:

        # if we've never seen this header - add it to the lookup dictionary
        # as an entry in a list. This means multiple entries can have the same header (which
        # we're saying is OK)
        if entry[0] not in lookup:
            lookup[entry[0]] = [entry[1]]
            updated.append(entry)

        # however, if we HAVE seen this header before..
        else:
            found_dupe = False
            
            # for each sequence associated wth that header
            # ask if we found a duplicate 
            for d in lookup[entry[0]]:
                if d == entry[1]:
                    found_dupe=True

            # if we found no duplicates add
            if not found_dupe:
                lookup[entry[0]].append(entry[1])
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
    seq_to_header: dict[str, str] = {}
    for entry in dataset:
        if entry[1] in seq_to_header:
            raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (seq_to_header[entry[1]], entry[0]))
        seq_to_header[entry[1]] = entry[0]



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
    lookup: set[str] = set()
    updated: list[list[str]] = []

    for entry in dataset:

        # if the sequence is in the lookup
        # set
        if entry[1] in lookup:
            continue

        # else its not in the lookup set so 
        # add it and then add to the updated list
        else:
            lookup.add(entry[1])
            updated.append(entry)
    return updated
            

    
                

                
