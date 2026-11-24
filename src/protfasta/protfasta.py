"""
protfasta - A simple but robust FASTA parser explicitly for protein sequences.

This module contains the internal workflow functions that orchestrate
duplicate handling and invalid-sequence processing for :func:`read_fasta`.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""

from __future__ import annotations

from typing import Optional, Union

from . import utilities as _utilities
from . import io as _io
from .protfasta_exceptions import ProtfastaException

    

####################################################################################################
#
#
def _deal_with_invalid_sequences(
    raw: list[list[str]],
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    verbose: bool = False,
    correction_dictionary: Optional[Union[dict[str, str], bool]] = None,
) -> list[list[str]]:
    """Process sequences that contain non-standard amino-acid characters.

    Depending on *invalid_sequence_action*, invalid residues are either
    left alone, cause an exception, are removed (whole sequence), or are
    converted to standard residues.

    Parameters
    ----------
    raw : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    invalid_sequence_action : str, optional
        Strategy for handling invalid residues.  One of:

        * ``'fail'`` (default) -- raise on the first invalid residue.
        * ``'remove'`` -- discard any sequence containing invalid
          residues.
        * ``'convert'`` -- convert using *correction_dictionary* (or
          built-in defaults), then fail if unconvertible residues
          remain.
        * ``'convert-ignore'`` -- convert what can be converted,
          silently keep any remaining invalid residues.
        * ``'ignore'`` -- accept invalid residues without changes.

    alignment : bool, optional
        When ``True``, dash characters (``'-'``) are treated as valid
        gap characters rather than invalid residues.  Default ``False``.

    verbose : bool, optional
        If ``True``, print progress information to stdout.
        Default ``False``.

    correction_dictionary : dict, bool, or None, optional
        Mapping of non-standard characters to their replacements.
        ``None`` uses the built-in table; ``False`` disables conversion
        (used internally during ``'convert-remove'`` two-pass logic).

    Returns
    -------
    list[list[str]]
        The (possibly modified) list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        If *invalid_sequence_action* is ``'fail'`` or ``'convert'`` and
        invalid residues are found (or remain after conversion).
    """

    # fail on the first invalid sequence
    if invalid_sequence_action == 'fail':
        _utilities.fail_on_invalid_sequences(raw, alignment)
        return raw

    # simply remove invalid sequences
    if invalid_sequence_action == 'remove':
        updated = _utilities.remove_invalid_sequences(raw, alignment)
        
        if verbose:
            print('[INFO]: Removed %i of %i due to sequences with invalid characters' % (len(raw) - len(updated), len(raw)))

        return updated

    # convert invalid sequences
    if invalid_sequence_action == 'convert' or invalid_sequence_action == 'convert-ignore':

        # correction_dictionary can arrive as the sentinel False (see the
        # 'convert-remove' two-pass logic in read_fasta); normalize it to None
        # here, which selects the built-in conversion table
        cd = correction_dictionary if isinstance(correction_dictionary, dict) else None
        (updated, count) = _utilities.convert_invalid_sequences(raw, cd, alignment)
        if verbose:
            print('[INFO]: Converted invalid residues in %i sequences' % (count))

        # note we then rescan in case there were still characters we couldn't deal with
        if invalid_sequence_action == 'convert':
            try:
                _utilities.fail_on_invalid_sequences(updated, alignment)
            except ProtfastaException as e:
                raise ProtfastaException("\n\n******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s"%(str(e)))

        return updated

    # ignore invalid sequences
    if invalid_sequence_action == 'ignore':
        return raw
        
    raise ProtfastaException("Invalid option passed to the selector 'invalid_sequence_action': %s" %(invalid_sequence_action))



####################################################################################################
#
#
def _deal_with_duplicate_records(
    raw: list[list[str]],
    duplicate_record_action: str = 'ignore',
    verbose: bool = False,
) -> list[list[str]]:
    """Process duplicate FASTA records (same header **and** sequence).

    Parameters
    ----------
    raw : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    duplicate_record_action : str, optional
        Strategy for handling duplicates.  One of:

        * ``'ignore'`` (default) -- keep all occurrences.
        * ``'fail'`` -- raise on the first duplicate.
        * ``'remove'`` -- keep only the first occurrence of each
          duplicate.

    verbose : bool, optional
        If ``True``, print the number of removed records to stdout.
        Default ``False``.

    Returns
    -------
    list[list[str]]
        The (possibly filtered) list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        If *duplicate_record_action* is ``'fail'`` and a duplicate
        record is found.
    """

    if duplicate_record_action == 'ignore':
        pass

    if duplicate_record_action == 'fail':
        _utilities.fail_on_duplicates(raw)

    if duplicate_record_action == 'remove':
        updated = _utilities.remove_duplicates(raw)
        if verbose:
            print('[INFO]: Removed %i of %i due to duplicate records ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw



####################################################################################################
#
#
def _deal_with_duplicate_sequences(
    raw: list[list[str]],
    duplicate_sequence_action: str = 'ignore',
    verbose: bool = False,
) -> list[list[str]]:
    """Process entries that share the same sequence (regardless of header).

    Parameters
    ----------
    raw : list[list[str]]
        Parsed FASTA data -- a list of ``[header, sequence]`` pairs.

    duplicate_sequence_action : str, optional
        Strategy for handling duplicate sequences.  One of:

        * ``'ignore'`` (default) -- keep all occurrences.
        * ``'fail'`` -- raise on the first duplicate sequence.
        * ``'remove'`` -- keep only the first occurrence.

    verbose : bool, optional
        If ``True``, print the number of removed entries to stdout.
        Default ``False``.

    Returns
    -------
    list[list[str]]
        The (possibly filtered) list of ``[header, sequence]`` pairs.

    Raises
    ------
    ProtfastaException
        If *duplicate_sequence_action* is ``'fail'`` and a duplicate
        sequence is found.
    """

    if duplicate_sequence_action == 'ignore':
        pass

    if duplicate_sequence_action == 'fail':
        _utilities.fail_on_duplicate_sequences(raw)
        
    if duplicate_sequence_action == 'remove':
        updated = _utilities.remove_duplicate_sequences(raw)
        if verbose:
            print('[INFO]: Removed %i of %i due to duplicate sequences ' % (len(raw) - len(updated), len(raw)))
        return updated

    return raw


####################################################################################################
#
#
def _deal_with_sequences_emptied_by_conversion(
    raw: list[list[str]],
    empty_sequence_action: str = 'fail',
    verbose: bool = False,
) -> list[list[str]]:
    """Apply *empty_sequence_action* to sequences that conversion left empty.

    The conversion table deletes some characters outright (``'*'``,
    ``'-'``, ``' '``), so a sequence made up of nothing else comes out of
    the invalid-residue step as an empty string.  That is the same thing
    as a header with no sequence, so it gets the same treatment: an error
    by default, or dropped/kept if the caller asked for that.

    Parameters
    ----------
    raw : list[list[str]]
        ``[header, sequence]`` pairs, after residue conversion.

    empty_sequence_action : str, optional
        ``'fail'`` (default) raises on the first empty sequence,
        ``'remove'`` drops such records and ``'ignore'`` keeps them.

    verbose : bool, optional
        If ``True``, print how many records were removed.  Default
        ``False``.

    Returns
    -------
    list[list[str]]
        The records that survive.  The input list is returned as it is
        when nothing had to be dropped.

    Raises
    ------
    ProtfastaException
        If *empty_sequence_action* is ``'fail'`` and a sequence is empty.
    """
    # the common case: nothing was emptied, so nothing needs rebuilding
    if all([entry[1] for entry in raw]):
        return raw

    updated = []
    for entry in raw:
        if entry[1] or _io._keep_empty_record(entry[0], empty_sequence_action, emptied_by_conversion=True):
            updated.append(entry)

    if verbose and empty_sequence_action == 'remove':
        print('[INFO]: Removed %i of %i records left with no sequence after conversion' % (len(raw) - len(updated), len(raw)))

    return updated
