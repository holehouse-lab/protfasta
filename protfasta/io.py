"""
protfasta - A simple but robust FASTA parser explicitly for protein sequences.

This module handles FASTA file I/O and validation of input arguments passed
to the public ``read_fasta`` API.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""

from __future__ import annotations

from typing import Callable, Optional

from .protfasta_exceptions import ProtfastaException


def check_inputs(
    expect_unique_header: bool,
    header_parser: Optional[Callable[[str], str]],
    check_header_parser: bool,
    duplicate_record_action: str,
    duplicate_sequence_action: str,
    invalid_sequence_action: str,
    alignment: bool,
    return_list: bool,
    output_filename: Optional[str],
    verbose: bool,
    correction_dictionary: Optional[dict[str, str]],
) -> None:
    """Validate all input arguments passed to :func:`read_fasta`.

    This is a stateless guard function: it either returns ``None`` when
    every argument is acceptable, or raises a
    :class:`~protfasta.protfasta_exceptions.ProtfastaException` describing
    the first problem it encounters.

    If new functionality is added to ``read_fasta``, the corresponding
    keyword must be validated here.

    Parameters
    ----------
    expect_unique_header : bool
        Whether every FASTA header in the file is expected to be unique.

    header_parser : callable or None
        A user-supplied function that accepts a single header string and
        returns a transformed header string.  ``None`` means no custom
        parsing.

    check_header_parser : bool
        When ``True``, ``header_parser`` is tested with a dummy string
        before the file is read to catch obvious problems early.

    duplicate_record_action : str
        How to handle duplicate records.  Must be one of
        ``'ignore'``, ``'fail'``, or ``'remove'``.

    duplicate_sequence_action : str
        How to handle duplicate sequences.  Must be one of
        ``'ignore'``, ``'fail'``, or ``'remove'``.

    invalid_sequence_action : str
        How to handle invalid amino-acid characters.  Must be one of
        ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``,
        ``'convert-ignore'``, or ``'convert-remove'``.

    alignment : bool
        Whether the file should be treated as an alignment (dashes are
        kept as valid gap characters).

    return_list : bool
        Whether the caller expects a list (``True``) or dict (``False``)
        return type.

    output_filename : str or None
        Optional path to write the final processed sequences to.

    verbose : bool
        Whether to emit informational messages to stdout.

    correction_dictionary : dict or None
        Optional mapping of non-standard characters to their
        replacements (e.g. ``{'B': 'N'}``).  Overrides the built-in
        conversion table when provided.

    Raises
    ------
    ProtfastaException
        If any argument fails validation.
    """

    # check the expect_unique_header keyword
    if type(expect_unique_header) != bool:
        raise ProtfastaException("keyword 'expect_unique_header' must be a boolean")

    # validate the header_parser
    if check_header_parser:
        if header_parser is not None:
            if not callable(header_parser):
                raise ProtfastaException("keyword 'header_parser' must be a function [tested with callable()]")
            
            try:
                tst_string = 'this test string should work'
                a = header_parser(tst_string)
                if type(a) != str:
                    raise ProtfastaException('Something went wrong when testing the header_parser function.\nFunction completed but return value was not a string')
            except Exception as e:
                raise ProtfastaException(f'Something went wrong when testing the header_parser function using string: {tst_string}.\nMaybe you should set check_header_parser to False? \nException: {e}')
            

    # check the duplicates_record_action 
    if duplicate_record_action not in ['ignore','fail','remove']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail','remove'")

    # check the duplicates_record_action 
    if duplicate_sequence_action not in ['ignore','fail','remove']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail', 'remove'")


    # check the invalid_sequence 
    if invalid_sequence_action not in ['ignore','fail','remove','convert','convert-ignore', 'convert-remove']:
        raise ProtfastaException("keyword 'invalid_sequence' must be one of 'ignore','fail','remove','convert','convert-ignore', 'convert-remove'")

    # check the return_list
    if type(return_list) != bool:
        raise ProtfastaException("keyword 'verbose' must be a boolean")

    # check the return_list
    if output_filename is not None:
        if type(output_filename) != str:
            raise ProtfastaException("keyword 'output_filename' must be a string")

    # check verbose
    if type(verbose) != bool:
        raise ProtfastaException("keyword 'verbose' must be a boolean")

    if type(alignment) != bool:
        raise ProtfastaException("keyword 'alignment' must be a boolean")
        
    if duplicate_record_action == 'ignore':
        if expect_unique_header is True:
            raise ProtfastaException('Cannot expect unique headers and ignore duplicate records')
            
    # checks correction dictionary
    if correction_dictionary is not None:
        if type(correction_dictionary) != dict:
            raise ProtfastaException("If provided, keyword 'correction_dictionary' must be a dictionary")





####################################################################################################
#
#
def internal_parse_fasta_file(
    filename: str,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
) -> list[list[str]]:
    """Low-level FASTA file parser.

    Reads a FASTA file from disk and returns its contents as a list of
    ``[header, sequence]`` pairs.  Header lines must begin with ``">"``
    and occupy a single line; sequence data may span multiple lines.

    This is an internal helper -- most callers should use
    :func:`protfasta.read_fasta` instead.

    Parameters
    ----------
    filename : str
        Absolute or relative path to a FASTA file.

    expect_unique_header : bool, optional
        If ``True`` (the default), a
        :class:`~protfasta.protfasta_exceptions.ProtfastaException` is
        raised when a duplicate header is encountered.

    header_parser : callable or None, optional
        A function that accepts a raw header string and returns a
        (possibly transformed) header string.  ``None`` means headers
        are used as-is (minus the leading ``">"``).  The function must
        accept exactly one ``str`` argument and return a ``str``.

    verbose : bool, optional
        If ``True``, informational messages are printed to stdout
        during parsing.

    Returns
    -------
    list[list[str]]
        A list of two-element lists ``[header, sequence]`` in the order
        they appear in the file.  Sequences are upper-cased.

    Raises
    ------
    ProtfastaException
        If the file cannot be found or a duplicate header is detected
        (when *expect_unique_header* is ``True``).
    """
        
    # read in the file...
    try:
        with open(filename,'r') as fh:
            content = fh.readlines()
    except FileNotFoundError:
        raise ProtfastaException('Unable to find file: %s'%(filename))
    
    if verbose:
        print('[INFO]: Read in file with %i lines'%(len(content)))

    # note, we'll keep the ability to directly parse dictionaries
    return _parse_fasta_all(content, expect_unique_header=expect_unique_header, header_parser=header_parser, verbose=verbose)
    


####################################################################################################
#
#
def _parse_fasta_all(
    content: list[str],
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
) -> list[list[str]]:
    """Parse raw FASTA-file lines into ``[header, sequence]`` pairs.

    This is the core parsing engine used by
    :func:`internal_parse_fasta_file`.  It operates on an already-read
    list of lines rather than a file path, making it easy to test in
    isolation.

    Parameters
    ----------
    content : list[str]
        Lines read from a FASTA file (e.g. via ``file.readlines()``).
        Each element is expected to be a single line, optionally
        terminated by a newline character.

    expect_unique_header : bool, optional
        If ``True`` (the default), a
        :class:`~protfasta.protfasta_exceptions.ProtfastaException` is
        raised when a duplicate header is encountered.

    header_parser : callable or None, optional
        A function ``(str) -> str`` used to transform raw header strings.
        ``None`` means headers are used verbatim (minus the leading
        ``">"`` character).

    verbose : bool, optional
        If ``True``, prints the number of recovered sequences to stdout.

    Returns
    -------
    list[list[str]]
        A list of two-element lists ``[header, sequence]``.  Sequences
        are upper-cased and concatenated from any multi-line runs in the
        input.

    Raises
    ------
    ProtfastaException
        If *expect_unique_header* is ``True`` and a duplicate header is
        found.
    """

    
    return_data=[]

    # note - using a dictionary for all_headers to look up
    # unique headers is MANY MANY MANY orders of magnitude faster than doing
    # this with a list of a set
    all_headers={}
    
    ## START OF PARSING FUNCTION    
    seq=''
    header=''

    def update():
    
        if header in all_headers:
            if expect_unique_header:
                raise ProtfastaException('Found duplicate header (%s)' %(header))
        else:
            all_headers[header] = True

        return_data.append([header, seq.upper()])


    for line in content:
        
        sline = line.strip()

        # if empty line just skip...
        if len(sline) == 0:
            continue

        # if  first non-whitespace character is a '>'
        if sline[0] == '>':

            # get the current header line
            h = sline[1:]
            
            # if we'd previously had a sequence assigned, means we have just started a 
            # 'new' sequence
            if len(seq) > 0:
                
                # if we're in the processing of parsing a string this means we found a new
                # header so update the return_data with header and sequence
                update()

                                    
            # reset the header and sequence
            if header_parser:
                header = header_parser(h)
            else:
                header = h
            seq=''
            
        else:            
            # we're on a line that is neither empty nor started with
            # a '>' so it's treated as a sequence line
            seq = seq + sline

    # if we exit with a sequence in toe, there's one final sequence to
    # add...        
    if len(seq) > 0:
        update()

    if verbose:
        print('[INFO]: Parsed file to recover %i sequences' %(len(return_data)))


    return return_data



