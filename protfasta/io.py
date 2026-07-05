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

from typing import Callable, Iterator, Optional, Union

from . import utilities as _utilities
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
def _iter_fasta(
    filename: str,
    header_parser: Optional[Callable[[str], str]] = None,
):
    """Yield raw ``(header, sequence)`` pairs from a FASTA file, streaming.

    Low-level streaming parse loop used internally by
    :func:`_stream_fasta` (and therefore by
    :func:`protfasta.read_fasta_stream`).  It keeps peak memory to
    ``O(single record)`` by reading the file line-by-line and joining
    each record's sequence lines once.

    No duplicate detection, invalid-residue handling, or alignment-gap
    logic is performed -- each record is yielded exactly as parsed.  The
    public streaming entry point, :func:`protfasta.read_fasta_stream`,
    layers that sanitization on top.

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


####################################################################################################
#
#
def _write_stream_record(fh, header: str, seq: str, linelength: int = 60) -> None:
    """Write a single ``[header, sequence]`` record to an open file handle.

    Used by :func:`_stream_fasta` to tee sanitized records to
    *output_filename* as they are produced.  Formatting matches
    :func:`protfasta.write_fasta` with its default ``linelength`` of 60
    residues per line and a blank separator line between records.

    Parameters
    ----------
    fh : file object
        An open, writeable text file handle.

    header : str
        The record header (written without the leading ``">"``, which is
        added here).

    seq : str
        The amino-acid sequence.  Must be non-empty.

    linelength : int, optional
        Number of residues per output line.  Default ``60``.

    Raises
    ------
    ProtfastaException
        If *seq* is empty (mirrors :func:`protfasta.write_fasta`).
    """
    if len(seq) < 1:
        raise ProtfastaException('Seqence associated with [%s] is empty' % (header))

    fh.write('>')
    fh.write(header)
    fh.write('\n')
    for start in range(0, len(seq), linelength):
        fh.write(seq[start:start + linelength])
        fh.write('\n')
    # Blank separator line between records.
    fh.write('\n')


####################################################################################################
#
#
def _stream_fasta(
    filename: str,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    duplicate_sequence_action: str = 'ignore',
    duplicate_record_action: str = 'fail',
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    return_list: bool = False,
    output_filename: Optional[str] = None,
    correction_dictionary: Optional[dict[str, str]] = None,
    verbose: bool = False,
) -> Iterator[Union[tuple[str, str], list[str]]]:
    """Stream a FASTA file record-by-record with full sanitization.

    This is the streaming engine behind :func:`protfasta.read_fasta_stream`.
    It applies the same processing steps as :func:`protfasta.read_fasta`
    -- header-uniqueness checks, duplicate-record and duplicate-sequence
    handling, and invalid-residue handling -- but yields one record at a
    time instead of materializing the whole dataset.

    The processing order matches :func:`protfasta.read_fasta`:

    1. Header uniqueness (*expect_unique_header*).
    2. Duplicate records (*duplicate_record_action*).
    3. Duplicate sequences (*duplicate_sequence_action*).
    4. Invalid residues (*invalid_sequence_action*).
    5. Optional tee to *output_filename*.

    Peak memory is ``O(number of records)`` for the auxiliary
    duplicate/uniqueness bookkeeping (headers plus 16-byte sequence
    digests -- never whole sequences) and ``O(single record)`` for the
    sequence data itself.  When *expect_unique_header* is ``False`` and
    all duplicate actions are ``'ignore'``, the auxiliary bookkeeping is
    skipped entirely and memory is flat regardless of file size.

    Parameters
    ----------
    filename : str
        Path to the FASTA file to read.

    expect_unique_header : bool, optional
        If ``True`` (default), raise on the first duplicate header.

    header_parser : callable or None, optional
        Optional ``(str) -> str`` transform applied to every raw header.

    duplicate_sequence_action : str, optional
        One of ``'ignore'``, ``'fail'``, or ``'remove'``.  Default
        ``'ignore'``.

    duplicate_record_action : str, optional
        One of ``'ignore'``, ``'fail'``, or ``'remove'``.  Default
        ``'fail'``.

    invalid_sequence_action : str, optional
        One of ``'ignore'``, ``'fail'``, ``'remove'``, ``'convert'``,
        ``'convert-ignore'``, or ``'convert-remove'``.  Default
        ``'fail'``.

    alignment : bool, optional
        If ``True``, dashes (``'-'``) are treated as valid gap
        characters.  Default ``False``.

    return_list : bool, optional
        If ``True``, yield ``[header, sequence]`` lists; otherwise yield
        ``(header, sequence)`` tuples (default).

    output_filename : str or None, optional
        If provided, each sanitized record is written to this path as it
        is yielded.  The output file is only complete once the generator
        has been fully consumed.

    correction_dictionary : dict or None, optional
        Custom character-replacement mapping used by the ``'convert'``
        actions.  ``None`` uses the built-in table.

    verbose : bool, optional
        If ``True``, emit an opening message and, when the generator is
        exhausted, a summary of removed/converted counts.

    Yields
    ------
    tuple[str, str] or list[str]
        ``(header, sequence)`` pairs (or ``[header, sequence]`` lists when
        *return_list* is ``True``) in file order.  Sequences are
        upper-cased and sanitized.

    Raises
    ------
    ProtfastaException
        On a duplicate header/record/sequence (for the relevant ``'fail'``
        actions) or an invalid residue (for ``'fail'``/``'convert'``).
        Because parsing is lazy, these are raised mid-iteration, at the
        offending record.
    """

    # Only allocate bookkeeping structures for the actions that need them.
    seen_headers: Optional[set[str]] = set() if expect_unique_header else None

    # duplicate records: header -> set of sequence digests
    record_lookup: Optional[dict[str, set[bytes]]] = (
        {} if duplicate_record_action in ('fail', 'remove') else None
    )

    # duplicate sequences: sequence digest -> first header seen (the header
    # is retained so the 'fail' message can name both offenders).
    seq_lookup: Optional[dict[bytes, str]] = (
        {} if duplicate_sequence_action in ('fail', 'remove') else None
    )

    need_digest = record_lookup is not None or seq_lookup is not None

    n_read = 0
    n_yielded = 0
    n_dup_records_removed = 0
    n_dup_seqs_removed = 0
    n_invalid_removed = 0
    n_converted = 0

    # Large write buffer to minimise syscall overhead, matching write_fasta.
    out_fh = open(output_filename, 'w', buffering=1024 * 1024) if output_filename else None

    if verbose:
        print('[INFO]: Streaming file %s' % (filename))

    try:
        for header, seq in _iter_fasta(filename, header_parser=header_parser):
            n_read += 1

            # 1. header uniqueness
            if seen_headers is not None:
                if header in seen_headers:
                    raise ProtfastaException('Found duplicate header (%s)' % (header))
                seen_headers.add(header)

            digest = _utilities._seq_hash(seq) if need_digest else b''

            # 2. duplicate records (identical header AND sequence)
            if record_lookup is not None:
                seen = record_lookup.get(header)
                if seen is None:
                    record_lookup[header] = {digest}
                elif digest in seen:
                    if duplicate_record_action == 'fail':
                        raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (header, seq))
                    n_dup_records_removed += 1
                    continue
                else:
                    seen.add(digest)

            # 3. duplicate sequences (identical sequence, any header)
            if seq_lookup is not None:
                if digest in seq_lookup:
                    if duplicate_sequence_action == 'fail':
                        raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (seq_lookup[digest], header))
                    n_dup_seqs_removed += 1
                    continue
                seq_lookup[digest] = header

            # 4. invalid-residue handling (per record)
            if invalid_sequence_action == 'ignore':
                pass

            elif invalid_sequence_action == 'fail':
                (status, info) = _utilities.check_sequence_is_valid(seq, alignment)
                if status is not True:
                    raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (info, header, seq))

            elif invalid_sequence_action == 'remove':
                (status, _info) = _utilities.check_sequence_is_valid(seq, alignment)
                if status is not True:
                    n_invalid_removed += 1
                    continue

            elif invalid_sequence_action == 'convert':
                new_seq = _utilities.convert_to_valid(seq, correction_dictionary, alignment)
                if new_seq != seq:
                    n_converted += 1
                (status, info) = _utilities.check_sequence_is_valid(new_seq, alignment)
                if status is not True:
                    inner = 'Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (info, header, new_seq)
                    raise ProtfastaException("\n\n******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s" % (inner))
                seq = new_seq

            elif invalid_sequence_action == 'convert-ignore':
                new_seq = _utilities.convert_to_valid(seq, correction_dictionary, alignment)
                if new_seq != seq:
                    n_converted += 1
                seq = new_seq

            elif invalid_sequence_action == 'convert-remove':
                new_seq = _utilities.convert_to_valid(seq, correction_dictionary, alignment)
                if new_seq != seq:
                    n_converted += 1
                (status, _info) = _utilities.check_sequence_is_valid(new_seq, alignment)
                if status is not True:
                    n_invalid_removed += 1
                    continue
                seq = new_seq

            # 5. tee the sanitized record to disk (if requested)
            if out_fh is not None:
                _write_stream_record(out_fh, header, seq)

            n_yielded += 1
            if return_list:
                yield [header, seq]
            else:
                yield (header, seq)

        if verbose:
            if duplicate_record_action == 'remove':
                print('[INFO]: Removed %i of %i due to duplicate records ' % (n_dup_records_removed, n_read))
            if duplicate_sequence_action == 'remove':
                print('[INFO]: Removed %i of %i due to duplicate sequences ' % (n_dup_seqs_removed, n_read))
            if invalid_sequence_action in ('convert', 'convert-ignore', 'convert-remove'):
                print('[INFO]: Converted %i sequences to valid sequences' % (n_converted))
            if invalid_sequence_action in ('remove', 'convert-remove'):
                print('[INFO]: Removed %i of %i due to sequences with invalid characters' % (n_invalid_removed, n_read))
            print('[INFO]: Streamed %i of %i records from %s' % (n_yielded, n_read, filename))

    finally:
        if out_fh is not None:
            out_fh.close()


