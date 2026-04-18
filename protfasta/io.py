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
        
    # Stream the file line-by-line rather than materializing the whole
    # file with readlines().  This keeps peak memory to O(single record)
    # which is essential for multi-gigabyte FASTA files.
    try:
        fh = open(filename, 'r')
    except FileNotFoundError:
        raise ProtfastaException('Unable to find file: %s' % (filename))

    if verbose:
        print('[INFO]: Read in file %s (streaming)' % (filename))

    with fh:
        return _parse_fasta_all(
            fh,
            expect_unique_header=expect_unique_header,
            header_parser=header_parser,
            verbose=verbose,
        )


####################################################################################################
#
#
def _parse_fasta_all(
    content,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
) -> list[list[str]]:
    """Parse FASTA content into ``[header, sequence]`` pairs.

    This is the core parsing engine used by
    :func:`internal_parse_fasta_file`.  It accepts any iterable of
    lines -- typically an open file handle (streamed, preferred for
    large files) or a pre-loaded ``list[str]`` -- making it easy to
    test in isolation.

    Parameters
    ----------
    content : Iterable[str]
        Iterable yielding lines from a FASTA file.  Each element is
        expected to be a single line, optionally terminated by a
        newline character.  An open file handle is accepted and will
        be consumed lazily.

    expect_unique_header : bool, optional
        If ``True`` (the default), a
        :class:`~protfasta.protfasta_exceptions.ProtfastaException` is
        raised when a duplicate header is encountered.  When ``False``
        no header tracking is performed (saves memory on large files).

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

    return_data: list[list[str]] = []

    # Only allocate a header-tracking set when we actually need it.
    # Using a set (not a dict-with-True-values) saves memory; skipping
    # it entirely when expect_unique_header is False saves even more
    # on files with hundreds of millions of records.
    seen_headers: Optional[set[str]] = set() if expect_unique_header else None

    # Accumulate sequence lines into a list and join once per record;
    # this is O(n) per sequence instead of the O(n^2) that repeated
    # string concatenation can degenerate into.
    seq_parts: list[str] = []
    header = ''
    have_record = False

    for line in content:

        # rstrip() handles \n, \r, and any trailing whitespace in one
        # C-level pass.  This is cheaper than a full strip() and a
        # subsequent test, since most FASTA lines have no leading WS.
        line = line.rstrip()

        if not line:
            continue

        if line[0] == '>':
            # Flush the previous record, but only if we accumulated at
            # least one sequence line (matching the original behaviour
            # of silently skipping headers with no sequence).
            if have_record and seq_parts:
                if seen_headers is not None:
                    if header in seen_headers:
                        raise ProtfastaException('Found duplicate header (%s)' % (header))
                    seen_headers.add(header)
                return_data.append([header, ''.join(seq_parts).upper()])

            # Start the new record.
            h = line[1:]
            header = header_parser(h) if header_parser else h
            seq_parts = []
            have_record = True
        else:
            seq_parts.append(line)

    # Flush the final record if the file didn't end with a blank line.
    if have_record and seq_parts:
        if seen_headers is not None:
            if header in seen_headers:
                raise ProtfastaException('Found duplicate header (%s)' % (header))
            seen_headers.add(header)
        return_data.append([header, ''.join(seq_parts).upper()])

    if verbose:
        print('[INFO]: Parsed file to recover %i sequences' % (len(return_data)))

    return return_data



####################################################################################################
#
#
def iter_fasta(
    filename: str,
    header_parser: Optional[Callable[[str], str]] = None,
):
    """Yield ``(header, sequence)`` pairs from a FASTA file, streaming.

    This is a memory-efficient alternative to :func:`protfasta.read_fasta`
    designed for very large files (hundreds of millions of sequences)
    where holding the entire dataset in memory is not feasible.

    No duplicate detection, invalid-residue handling, or alignment-gap
    logic is performed -- each record is yielded as parsed.  Callers
    that need those features should consume the iterator and apply
    their own filtering.

    Parameters
    ----------
    filename : str
        Path to a FASTA file.

    header_parser : callable or None, optional
        Optional ``(str) -> str`` transform applied to every raw header.

    Yields
    ------
    tuple[str, str]
        ``(header, sequence)`` pairs in the order they appear in the
        file.  Sequences are upper-cased.

    Raises
    ------
    ProtfastaException
        If the file cannot be opened.
    """
    try:
        fh = open(filename, 'r')
    except FileNotFoundError:
        raise ProtfastaException('Unable to find file: %s' % (filename))

    seq_parts: list[str] = []
    header = ''
    have_record = False

    try:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue

            if line[0] == '>':
                if have_record and seq_parts:
                    yield (header, ''.join(seq_parts).upper())
                h = line[1:]
                header = header_parser(h) if header_parser else h
                seq_parts = []
                have_record = True
            else:
                seq_parts.append(line)

        if have_record and seq_parts:
            yield (header, ''.join(seq_parts).upper())
    finally:
        fh.close()


