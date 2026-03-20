"""
protfasta - A simple but robust FASTA parser explicitly for protein sequences.

This package provides two user-facing functions:

* :func:`read_fasta` -- read and sanitize a FASTA file.
* :func:`write_fasta` -- write sequences to a FASTA file.
"""

from __future__ import annotations

import os
from typing import Callable, Optional, Union

from protfasta import utilities as _utilities
from protfasta import io as _io
from protfasta._configs import STANDARD_AAS, STANDARD_CONVERSION
from protfasta import protfasta as _protfasta
from protfasta.protfasta_exceptions import ProtfastaException

## ------------------------------------------------------------
# READTHEDOCS versioning hack
#
# Generate _version.py if missing and in the Read the Docs environment
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile('../protfasta/_version.py'):   
    import versioningit            
    __version__ = versioningit.get_version('../')
else:
    from ._version import __version__
    
    
## ------------------------------------------------------------




_ROOT = os.path.abspath(os.path.dirname(__file__))


def _get_data(path: str) -> str:
    """Return the absolute path to a bundled data directory.

    Used internally to locate test data shipped with the package.

    Parameters
    ----------
    path : str
        Relative path (under the package ``data/`` directory) to resolve.

    Returns
    -------
    str
        Absolute filesystem path.
    """
    return os.path.join(_ROOT, 'data', path)


def read_fasta(
    filename: str,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    check_header_parser: bool = True,
    duplicate_sequence_action: str = 'ignore',
    duplicate_record_action: str = 'fail',
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    return_list: bool = False,
    output_filename: Optional[str] = None,
    correction_dictionary: Optional[dict[str, str]] = None,
    verbose: bool = False,
) -> Union[dict[str, str], list[list[str]]]:
    """Read a FASTA file, sanitize sequences, and return a dict or list.

    This is the primary entry point for **protfasta**.  At its simplest::

        sequences = read_fasta('proteins.fasta')

    returns a dictionary whose keys are FASTA headers and whose values
    are amino-acid sequences.  Many optional parameters allow automatic
    handling of duplicates, invalid residues, alignment gap characters,
    and more.

    Sanitization is applied in the following order:

    1. File is read, custom headers are parsed, and header uniqueness is
       checked (when *expect_unique_header* is ``True``).
    2. Duplicate records are processed (*duplicate_record_action*).
    3. Duplicate sequences are processed (*duplicate_sequence_action*).
    4. Invalid residues are processed (*invalid_sequence_action*).
    5. Final sequences are optionally written to *output_filename*.
    6. A dictionary or list is returned to the caller.

    Incompatible option combinations are caught before the file is read.

    Parameters
    ----------
    filename : str
        Path to the FASTA file to read.

    expect_unique_header : bool, optional
        If ``True`` (default), an exception is raised when a duplicate
        header is encountered during parsing.  Set to ``False`` when
        the file is known to contain duplicate headers -- in that case
        *return_list* should typically be ``True`` as well so that no
        entries are silently lost via dictionary-key overwriting.

    header_parser : callable or None, optional
        A function ``(str) -> str`` applied to every raw header before
        any uniqueness checks.  Useful for extracting accession IDs
        from structured headers.  When *check_header_parser* is
        ``True`` (the default) the function is smoke-tested with the
        string ``'this test string should work'`` before parsing
        begins.

    check_header_parser : bool, optional
        If ``True`` (default), *header_parser* is tested with a dummy
        string before the file is read to catch obvious problems early.
        Set to ``False`` to skip this pre-check.

    duplicate_record_action : str, optional
        How to handle records that are identical in both header **and**
        sequence.  Default ``'fail'``.

        * ``'ignore'``  -- keep all occurrences (requires
          *expect_unique_header* = ``False``).
        * ``'fail'``    -- raise an exception.
        * ``'remove'``  -- keep only the first occurrence.

    duplicate_sequence_action : str, optional
        How to handle entries that share the same sequence regardless of
        header.  Default ``'ignore'``.

        * ``'ignore'``  -- keep all occurrences.
        * ``'fail'``    -- raise an exception.
        * ``'remove'``  -- keep only the first occurrence.

    invalid_sequence_action : str, optional
        How to handle sequences containing non-standard amino-acid
        characters.  Default ``'fail'``.

        * ``'ignore'``         -- silently accept invalid residues.
        * ``'fail'``           -- raise an exception.
        * ``'remove'``         -- discard the entire sequence.
        * ``'convert'``        -- convert non-standard residues using
          *correction_dictionary* (or built-in defaults); fail if any
          unconvertible residues remain.
        * ``'convert-ignore'`` -- convert what can be converted, then
          ignore any remaining invalid residues.
        * ``'convert-remove'`` -- convert what can be converted, then
          discard sequences that still contain invalid residues.

    alignment : bool, optional
        If ``True``, dash (``'-'``) characters are treated as valid gap
        characters and are neither flagged as invalid nor converted.
        Default ``False``.

    return_list : bool, optional
        If ``True``, return a list of ``[header, sequence]`` pairs
        instead of a dictionary.  Required when duplicate headers are
        present and you want to keep all of them.  Default ``False``.

    output_filename : str or None, optional
        If provided, the final (sanitized) set of sequences is written
        to a new FASTA file at this path before the function returns.

    correction_dictionary : dict or None, optional
        A mapping of non-standard characters to replacement strings used
        when *invalid_sequence_action* involves conversion.  When
        ``None``, the built-in table is used:

        * ``B`` -> ``N``, ``U`` -> ``C``, ``X`` -> ``G``, ``Z`` -> ``Q``
        * ``*`` -> ``''``, ``-`` -> ``''``, ``' '`` -> ``''``

        A custom dictionary **replaces** the built-in table entirely.

    verbose : bool, optional
        If ``True``, informational messages are printed to stdout
        during each processing step.  Default ``False``.

    Returns
    -------
    dict[str, str] or list[list[str]]
        When *return_list* is ``False`` (default), a dictionary mapping
        headers to sequences.  When ``True``, a list of two-element
        lists ``[header, sequence]``.  Ordering always matches the
        original file.

    Raises
    ------
    ProtfastaException
        If any validation check fails or incompatible options are
        provided.
    """

    # first we sanity check all of the inputs provided. NOTE. If additional functionality is added, new
    # keywords MUST be sanity checked in this function
    _io.check_inputs(expect_unique_header,
                     header_parser, 
                     check_header_parser,
                     duplicate_record_action,
                     duplicate_sequence_action,
                     invalid_sequence_action, 
                     alignment,
                     return_list, 
                     output_filename,
                     verbose,
                     correction_dictionary)
    

    # the actual file i/o happens here
    raw = _io.internal_parse_fasta_file(filename, expect_unique_header=expect_unique_header, header_parser=header_parser, verbose=verbose)

    # first deal with duplicate records
    updated = _protfasta._deal_with_duplicate_records(raw, duplicate_record_action, verbose)

    # deal with duplicate sequences
    updated = _protfasta._deal_with_duplicate_sequences(updated, duplicate_sequence_action, verbose)

    # next decide how we deal with invalid amino acid sequences

    ##
    ## If we're using the convert-remove action...
    if invalid_sequence_action == 'convert-remove':

        # first run a convert ignore
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          'convert-ignore',
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=correction_dictionary)

        # then a remove on those that are left
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          'remove',
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=False)

    else:
        updated = _protfasta._deal_with_invalid_sequences(updated, 
                                                          invalid_sequence_action,
                                                          alignment=alignment,
                                                          verbose=verbose, 
                                                          correction_dictionary=correction_dictionary)
        

    # If we wanted to write the final set of sequences we're going to use...:
    if output_filename:
        write_fasta(updated, output_filename)

    # if we asked for a list...
    if return_list is True:
        pass    
    else:
        updated = _utilities.convert_list_to_dictionary(updated, verbose)

        
    return updated
        


# ------------------------------------------------------------------
#
def write_fasta(
    fasta_data: Union[dict[str, str], list[list[str]]],
    filename: str,
    linelength: Union[int, bool, None] = 60,
    verbose: bool = False,
    append_to_fasta: bool = False,
) -> None:
    """Write sequences to a FASTA file.

    Accepts sequence data as either a dictionary (header -> sequence) or
    a list of ``[header, sequence]`` pairs and writes a standards-
    compliant FASTA file to *filename*.

    Parameters
    ----------
    fasta_data : dict[str, str] or list[list[str]]
        Sequence data.  If a dictionary, keys are headers and values are
        amino-acid sequences.  If a list, each element must be a
        two-element list ``[header, sequence]``.

    filename : str
        Destination file path.  Should conventionally end with
        ``.fasta`` or ``.fa``, but this is not enforced.

    linelength : int, bool, or None, optional
        Maximum number of residues per line in the output.  Default is
        ``60`` (the UniProt convention).  Values below ``5`` are clamped
        to ``5``.  Set to ``0``, ``None``, or ``False`` to write each
        sequence on a single line.

    verbose : bool, optional
        Currently unused; reserved for future diagnostic output.
        Default ``False``.

    append_to_fasta : bool, optional
        If ``True``, new entries are appended to *filename* if it
        already exists; otherwise the file is created.  If ``False``
        (default), any existing file is overwritten.

    Returns
    -------
    None

    Raises
    ------
    ProtfastaException
        If a sequence is empty or a list element does not contain
        exactly two items.
    """

    # This part of code means we can pass either a dictionary or a list of lists
    # in for write_fasta to deal with

    if type(fasta_data) == list:
        def get_sequence():
            return (entry[0], entry[1])

        # quick validate
        for i in fasta_data:
            if len(i)  != 2:
                raise ProtfastaException('While processing a list for write_fasta_file at least one of the elements was not a 2-position sublist:\n%s'%(str(i)))

    if type(fasta_data) == dict:
        def get_sequence():
            return (entry, fasta_data[entry])
        

    # override line length for sane input. N
    if linelength == False or linelength == None or linelength < 1:
        linelength = False
        
    else:        
        # cast linelength to an integer here as a soft type checking, and
        # if it's shorter than 5 then reset to 5 
        linelength = int(linelength)
        if linelength < 5:
            linelength = 5

    # set the 'mode' for open. If append_to_file==False, use 'w' and overwrite
    # existing .fasta file. Otherwise use 'a' and add to existing file if it exists.
    if append_to_fasta==False:
        open_mode='w'
    else:
        open_mode='a'
    
    # open the file handle.
    with open(filename, open_mode) as fh:

        # for each entry
        for entry in fasta_data:

            (header, seq) = get_sequence()
            if len(seq) < 1:
                raise ProtfastaException('Seqence associated with [%s] is empty'%(header))
                

            # write the header line
            fh.write('>%s\n'%(header))
            
            # the $wrotenewline bool is ONLY here to
            # avoid the unlikely scenario in which the last character
            # is also an integer number of linelength, such that you'd
            # get TWO spaces between sequences. This ensures there is ALWAYS
            # only one blank line between sequences
            for i in range(0,len(seq)):
                fh.write('%s'%seq[i])
                wrotenewline=False

                # if linelength is valid
                if linelength:

                    # if we reach an integer number of line-length write
                    # a newline character
                    if (i+1) % linelength == 0:
                        fh.write('\n')
                        wrotenewline=True
            
            if wrotenewline:
                fh.write('\n')
            else:
                fh.write('\n\n')
    
