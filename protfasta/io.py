"""
protfasta - A simple but robust FASTA parser explicitly for protein sequences.

This module handles FASTA file I/O and validation of input arguments passed
to the public ``read_fasta`` / ``read_fasta_stream`` API.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""

from __future__ import annotations

import gc
import os
from typing import IO, Callable, Iterable, Iterator, Optional, Union

from . import utilities as _utilities
from .protfasta_exceptions import ProtfastaException


# ---------------------------------------------------------------------------
# Encoding policy
#
# FASTA files are read as UTF-8 and written as UTF-8, independent of the
# platform locale (which is what Python's open() would otherwise use, and
# which on Windows is typically not UTF-8). Two refinements:
#
#   * Reading uses 'utf-8-sig', which strips a leading byte-order mark if
#     one is present and is otherwise identical to 'utf-8'. Without this a
#     BOM (which Windows editors add routinely) hides the first '>' and
#     silently drops the first record.
#
#   * Both directions use the 'surrogateescape' error handler, so bytes
#     that are not valid UTF-8 (a Latin-1 accent in a header, say) are
#     carried through as lone surrogate code points rather than raising a
#     UnicodeDecodeError, and are written back out as the original bytes.
#     A stray byte inside a *sequence* therefore surfaces as an invalid
#     residue and is handled by invalid_sequence_action like any other.
# ---------------------------------------------------------------------------
_READ_ENCODING = 'utf-8-sig'
_WRITE_ENCODING = 'utf-8'
_ENCODING_ERRORS = 'surrogateescape'

# Large write buffer (1 MiB) to minimise syscall overhead on big outputs.
_WRITE_BUFFER = 1024 * 1024

# Type alias for anything we accept as a path.
PathLike = Union[str, 'os.PathLike[str]']


####################################################################################################
#
#
def check_filename(filename: object) -> None:
    """Validate that *filename* is something we can safely open.

    ``open()`` happily accepts an integer and interprets it as an already-open
    file descriptor, so a mistyped call such as ``read_fasta(0)`` would
    silently read from stdin rather than reporting a bad argument.  This guard
    restricts the input to a path -- a string or any
    :class:`os.PathLike` (e.g. :class:`pathlib.Path`).

    Parameters
    ----------
    filename : object
        The value passed as a FASTA filename.

    Raises
    ------
    ProtfastaException
        If *filename* is not a string or path-like object.
    """
    if not isinstance(filename, (str, os.PathLike)):
        raise ProtfastaException("keyword 'filename' must be a string or path-like object (got %s)" % (type(filename).__name__))


####################################################################################################
#
#
def _open_fasta(filename: PathLike) -> IO[str]:
    """Open a FASTA file for reading, converting OS errors into protfasta errors.

    Every failure mode of :func:`open` -- a missing file, a directory passed
    in place of a file, a permissions problem -- is surfaced as a
    :class:`~protfasta.protfasta_exceptions.ProtfastaException` so that
    callers only ever need to catch one exception type.  See the encoding
    policy at the top of this module for how the file is decoded.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the file to open.

    Returns
    -------
    file object
        An open, readable text file handle.

    Raises
    ------
    ProtfastaException
        If the file cannot be opened for any reason.
    """
    try:
        return open(filename, 'r', encoding=_READ_ENCODING, errors=_ENCODING_ERRORS)
    except FileNotFoundError:
        raise ProtfastaException('Unable to find file: %s' % (filename))
    except OSError as e:
        raise ProtfastaException('Unable to read file: %s\nException: %s' % (filename, e))


####################################################################################################
#
#
def _open_output(filename: PathLike, append: bool = False) -> IO[str]:
    """Open a FASTA file for writing, converting OS errors into protfasta errors.

    Parameters
    ----------
    filename : str or os.PathLike
        Destination path.

    append : bool, optional
        If ``True`` open in append mode, otherwise truncate/create.
        Default ``False``.

    Returns
    -------
    file object
        An open, writeable text file handle with a 1 MiB buffer.

    Raises
    ------
    ProtfastaException
        If the file cannot be opened (missing directory, permissions, ...).
    """
    try:
        return open(filename, 'a' if append else 'w', buffering=_WRITE_BUFFER, encoding=_WRITE_ENCODING, errors=_ENCODING_ERRORS)
    except OSError as e:
        raise ProtfastaException('Unable to open file for writing: %s\nException: %s' % (filename, e))


####################################################################################################
#
#
def check_inputs(
    expect_unique_header: bool,
    header_parser: Optional[Callable[[str], str]],
    check_header_parser: bool,
    duplicate_record_action: str,
    duplicate_sequence_action: str,
    invalid_sequence_action: str,
    alignment: bool,
    return_list: bool,
    output_filename: Optional[PathLike],
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
        parsing.  Must be callable whenever it is not ``None``.

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

    output_filename : str, os.PathLike, or None
        Optional path to write the final processed sequences to.

    verbose : bool
        Whether to emit informational messages to stdout.

    correction_dictionary : dict or None
        Optional mapping of non-standard characters to their
        replacements (e.g. ``{'B': 'N'}``).  Overrides the built-in
        conversion table when provided.  Keys must be non-empty strings
        and values must be strings.

    Raises
    ------
    ProtfastaException
        If any argument fails validation.
    """

    # check the expect_unique_header keyword
    if not isinstance(expect_unique_header, bool):
        raise ProtfastaException("keyword 'expect_unique_header' must be a boolean")

    if not isinstance(check_header_parser, bool):
        raise ProtfastaException("keyword 'check_header_parser' must be a boolean")

    # validate the header_parser. Note the callable() check applies whether
    # or not check_header_parser is set - that flag only controls the smoke
    # test with a dummy string, and a non-callable can never work
    if header_parser is not None:
        if not callable(header_parser):
            raise ProtfastaException("keyword 'header_parser' must be a function [tested with callable()]")

        if check_header_parser:
            tst_string = 'this test string should work'
            try:
                a = header_parser(tst_string)
            except Exception as e:
                raise ProtfastaException(f'Something went wrong when testing the header_parser function using string: {tst_string}.\nMaybe you should set check_header_parser to False? \nException: {e}')

            if not isinstance(a, str):
                raise ProtfastaException('Something went wrong when testing the header_parser function.\nFunction completed but return value was not a string (got %s)' % (type(a).__name__))

    # check the duplicate_record_action
    if duplicate_record_action not in ['ignore', 'fail', 'remove']:
        raise ProtfastaException("keyword 'duplicate_record_action' must be one of 'ignore','fail','remove'")

    # check the duplicate_sequence_action
    if duplicate_sequence_action not in ['ignore', 'fail', 'remove']:
        raise ProtfastaException("keyword 'duplicate_sequence_action' must be one of 'ignore','fail', 'remove'")

    # check the invalid_sequence_action
    if invalid_sequence_action not in ['ignore', 'fail', 'remove', 'convert', 'convert-ignore', 'convert-remove']:
        raise ProtfastaException("keyword 'invalid_sequence_action' must be one of 'ignore','fail','remove','convert','convert-ignore', 'convert-remove'")

    # check the return_list
    if not isinstance(return_list, bool):
        raise ProtfastaException("keyword 'return_list' must be a boolean")

    # check the output_filename
    if output_filename is not None:
        if not isinstance(output_filename, (str, os.PathLike)):
            raise ProtfastaException("keyword 'output_filename' must be a string or path-like object")

    # check verbose
    if not isinstance(verbose, bool):
        raise ProtfastaException("keyword 'verbose' must be a boolean")

    if not isinstance(alignment, bool):
        raise ProtfastaException("keyword 'alignment' must be a boolean")

    if duplicate_record_action == 'ignore':
        if expect_unique_header is True:
            raise ProtfastaException('Cannot expect unique headers and ignore duplicate records')

    # checks correction dictionary
    if correction_dictionary is not None:
        if not isinstance(correction_dictionary, dict):
            raise ProtfastaException("If provided, keyword 'correction_dictionary' must be a dictionary")
        _utilities._validate_correction_dictionary(correction_dictionary)


####################################################################################################
#
#
def _iter_fasta_lines(
    content: Iterable[str],
    header_parser: Optional[Callable[[str], str]] = None,
) -> Iterator[tuple[str, str]]:
    """Core FASTA parsing engine: turn an iterable of lines into records.

    This is the single parser behind everything in protfasta -- both the
    load-it-all path (:func:`_parse_fasta_all`) and the streaming path
    (:func:`_stream_fasta`) consume it, so the two can never disagree on
    what constitutes a record.

    It accepts any iterable of lines -- typically an open text-mode file
    handle (consumed lazily, so peak memory is one record) or a
    pre-loaded ``list[str]`` for testing.  Each element is expected to
    be one line, optionally newline-terminated.  Header lines begin with
    ``">"`` and occupy a single line; sequence data may span multiple
    lines.  Trailing whitespace is stripped from every line, blank lines
    are skipped, and anything before the first header is ignored.

    A header that is not followed by any sequence data is silently
    skipped, matching the behaviour protfasta has always had.

    Parameters
    ----------
    content : Iterable[str]
        Iterable yielding the lines of a FASTA file.

    header_parser : callable or None, optional
        A function ``(str) -> str`` applied to every raw header (with the
        leading ``">"`` already removed).  ``None`` means headers are used
        verbatim.

    Yields
    ------
    tuple[str, str]
        ``(header, sequence)`` pairs in file order.  Sequences are
        upper-cased and concatenated from any multi-line runs.

    Raises
    ------
    ProtfastaException
        If *header_parser* raises on a header encountered in the file.
    """

    # Accumulate sequence lines into a list and join once per record;
    # this is O(n) per sequence instead of the O(n^2) that repeated
    # string concatenation can degenerate into.
    seq_parts: list[str] = []
    header: Optional[str] = None

    for line in content:

        # rstrip() handles \n, \r, and any trailing whitespace in one
        # C-level pass.  This is cheaper than a full strip() and a
        # subsequent test, since most FASTA lines have no leading WS.
        line = line.rstrip()

        if not line:
            continue

        if line[0] == '>':
            # Flush the previous record, but only if we accumulated at
            # least one sequence line.
            if header is not None and seq_parts:
                yield (header, ''.join(seq_parts).upper())

            # Start the new record.
            h = line[1:]
            if header_parser is not None:
                try:
                    h = header_parser(h)
                except Exception as e:
                    raise ProtfastaException('header_parser raised an exception on header [%s]\nException: %s' % (_utilities._printable(h), e))
            header = h
            seq_parts = []
        else:
            seq_parts.append(line)

    # Flush the final record.
    if header is not None and seq_parts:
        yield (header, ''.join(seq_parts).upper())


####################################################################################################
#
#
def _parse_fasta_all(
    content: Iterable[str],
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
) -> list[list[str]]:
    """Parse FASTA content into a list of ``[header, sequence]`` pairs.

    A thin wrapper over :func:`_iter_fasta_lines` that materialises every
    record and, when asked, checks header uniqueness as it goes.  Used by
    :func:`internal_parse_fasta_file`; also handy for tests since it
    accepts a plain list of lines.

    Parameters
    ----------
    content : Iterable[str]
        Iterable yielding lines from a FASTA file (see
        :func:`_iter_fasta_lines`).

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
        A list of two-element lists ``[header, sequence]``.

    Raises
    ------
    ProtfastaException
        If *expect_unique_header* is ``True`` and a duplicate header is
        found.
    """

    return_data: list[list[str]] = []

    # Only allocate a header-tracking set when we actually need it; the set
    # holds references to the header strings already in return_data, so it
    # costs one pointer-sized slot per record rather than a copy.
    seen_headers: Optional[set[str]] = set() if expect_unique_header else None

    # Every [header, sequence] pair is a GC-tracked container, so a large
    # file makes the cyclic garbage collector run thousands of times to scan
    # objects that cannot possibly be part of a cycle - about 20% of the
    # parse time at ten million records. The loop below is bounded and
    # creates no cycles, so the collector is paused for its duration and
    # restored (to whatever state it was in) on the way out, exceptions
    # included. The streaming parser does not do this, because user code
    # runs between its yields.
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        for header, seq in _iter_fasta_lines(content, header_parser):
            if seen_headers is not None:
                if header in seen_headers:
                    raise ProtfastaException('Found duplicate header (%s)' % (_utilities._printable(header)))
                seen_headers.add(header)
            return_data.append([header, seq])
    finally:
        if gc_was_enabled:
            gc.enable()

    if verbose:
        print('[INFO]: Parsed file to recover %i sequences' % (len(return_data)))

    return return_data


####################################################################################################
#
#
def internal_parse_fasta_file(
    filename: PathLike,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
) -> list[list[str]]:
    """Low-level FASTA file parser.

    Reads a FASTA file from disk and returns its contents as a list of
    ``[header, sequence]`` pairs.  The file is streamed line-by-line
    rather than read into memory in one go, so peak memory is the parsed
    result plus one record, not the parsed result plus the raw file.

    This is an internal helper -- most callers should use
    :func:`protfasta.read_fasta` instead.

    Parameters
    ----------
    filename : str or os.PathLike
        Absolute or relative path to a FASTA file.

    expect_unique_header : bool, optional
        If ``True`` (the default), a
        :class:`~protfasta.protfasta_exceptions.ProtfastaException` is
        raised when a duplicate header is encountered.

    header_parser : callable or None, optional
        A function that accepts a raw header string and returns a
        (possibly transformed) header string.  ``None`` means headers
        are used as-is (minus the leading ``">"``).

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
        If the file cannot be opened or a duplicate header is detected
        (when *expect_unique_header* is ``True``).
    """

    fh = _open_fasta(filename)

    if verbose:
        print('[INFO]: Read in file %s (streaming)' % (filename))

    with fh:
        return _parse_fasta_all(fh,
                                expect_unique_header=expect_unique_header,
                                header_parser=header_parser,
                                verbose=verbose)


####################################################################################################
#
#
def _iter_fasta(
    filename: PathLike,
    header_parser: Optional[Callable[[str], str]] = None,
) -> Iterator[tuple[str, str]]:
    """Yield raw ``(header, sequence)`` pairs from a FASTA file, streaming.

    A convenience wrapper that opens *filename* and hands it to
    :func:`_iter_fasta_lines`.  No duplicate detection, invalid-residue
    handling, or alignment-gap logic is performed -- each record is
    yielded exactly as parsed.  The public streaming entry point,
    :func:`protfasta.read_fasta_stream`, layers that sanitization on top.

    Parameters
    ----------
    filename : str or os.PathLike
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
    with _open_fasta(filename) as fh:
        yield from _iter_fasta_lines(fh, header_parser)


####################################################################################################
#
#
def _format_record(header: str, seq: str, linelength: int = 60) -> str:
    """Render one FASTA record as text.

    Shared by :func:`protfasta.write_fasta` and the streaming tee in
    :func:`_stream_fasta`, so the two always produce byte-identical
    output.  The sequence is sliced into *linelength* chunks and the
    whole record (header line, sequence lines, blank separator line) is
    assembled with a single ``str.join`` so the caller can emit it with
    one ``write`` call.

    Parameters
    ----------
    header : str
        The record header (the leading ``">"`` is added here).

    seq : str
        The amino-acid sequence.

    linelength : int, optional
        Number of residues per output line.  Any falsy value (``0``,
        ``None``, ``False``) writes the sequence on a single line.
        Default ``60``.

    Returns
    -------
    str
        The formatted record, ending in a blank line.
    """
    if linelength:
        body = '\n'.join([seq[start:start + linelength] for start in range(0, len(seq), linelength)])
    else:
        body = seq
    return '>%s\n%s\n\n' % (header, body)


####################################################################################################
#
#
def _stream_fasta(
    filename: PathLike,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    duplicate_sequence_action: str = 'ignore',
    duplicate_record_action: str = 'fail',
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    return_list: bool = False,
    output_filename: Optional[PathLike] = None,
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
    duplicate/uniqueness bookkeeping (headers plus 16-byte digests --
    never whole sequences) and ``O(single record)`` for the sequence
    data itself.  When *expect_unique_header* is ``False`` and all
    duplicate actions are ``'ignore'``, the auxiliary bookkeeping is
    skipped entirely and memory is flat regardless of file size.

    The input file is opened before the output file, so a missing input
    never leaves an empty output file behind.

    Parameters
    ----------
    filename : str or os.PathLike
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

    output_filename : str, os.PathLike, or None, optional
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

    # duplicate records: a set of 16-byte digests of (header, sequence).
    # When headers are required to be unique a duplicate record (same
    # header AND sequence) is impossible - the header check fires first -
    # so the record check would only ever hash every sequence for nothing.
    seen_records: Optional[set[bytes]] = (
        set() if (duplicate_record_action in ('fail', 'remove') and not expect_unique_header) else None
    )

    # duplicate sequences: sequence digest -> first header seen (the header
    # is retained so the 'fail' message can name both offenders).
    seq_lookup: Optional[dict[bytes, str]] = (
        {} if duplicate_sequence_action in ('fail', 'remove') else None
    )

    # Build the converter once (it caches its translate table) rather than
    # once per record.
    convert = _utilities._make_converter(correction_dictionary, alignment) if invalid_sequence_action.startswith('convert') else None
    check_valid = _utilities.check_sequence_is_valid
    printable = _utilities._printable

    n_read = 0
    n_yielded = 0
    n_dup_records_removed = 0
    n_dup_seqs_removed = 0
    n_invalid_removed = 0
    n_converted = 0

    # Open the input first so a missing input file never creates an empty
    # output file.
    in_fh = _open_fasta(filename)
    out_fh: Optional[IO[str]] = None

    try:
        if output_filename is not None:
            out_fh = _open_output(output_filename)

        if verbose:
            print('[INFO]: Streaming file %s' % (filename))

        for header, seq in _iter_fasta_lines(in_fh, header_parser):
            n_read += 1

            # 1. header uniqueness
            if seen_headers is not None:
                if header in seen_headers:
                    raise ProtfastaException('Found duplicate header (%s)' % (printable(header)))
                seen_headers.add(header)

            # 2. duplicate records (identical header AND sequence)
            if seen_records is not None:
                key = _utilities._record_hash(header, seq)
                if key in seen_records:
                    if duplicate_record_action == 'fail':
                        raise ProtfastaException('Found duplicate entries of the following record\n:>%s\n%s' % (printable(header), printable(seq)))
                    n_dup_records_removed += 1
                    continue
                seen_records.add(key)

            # 3. duplicate sequences (identical sequence, any header)
            if seq_lookup is not None:
                digest = _utilities._seq_hash(seq)
                if digest in seq_lookup:
                    if duplicate_sequence_action == 'fail':
                        raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (printable(seq_lookup[digest]), printable(header)))
                    n_dup_seqs_removed += 1
                    continue
                seq_lookup[digest] = header

            # 4. invalid-residue handling (per record)
            if invalid_sequence_action == 'ignore':
                pass

            elif invalid_sequence_action == 'fail':
                (status, info) = check_valid(seq, alignment)
                if status is not True:
                    raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (printable(str(info)), printable(header), printable(seq)))

            elif invalid_sequence_action == 'remove':
                (status, _info) = check_valid(seq, alignment)
                if status is not True:
                    n_invalid_removed += 1
                    continue

            else:
                # one of the three convert-* actions
                assert convert is not None
                new_seq = convert(seq)
                if new_seq is not seq and new_seq != seq:
                    n_converted += 1
                seq = new_seq

                if invalid_sequence_action == 'convert':
                    (status, info) = check_valid(seq, alignment)
                    if status is not True:
                        inner = 'Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (printable(str(info)), printable(header), printable(seq))
                        raise ProtfastaException("\n\n******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s" % (inner))

                elif invalid_sequence_action == 'convert-remove':
                    (status, _info) = check_valid(seq, alignment)
                    if status is not True:
                        n_invalid_removed += 1
                        continue

                # 'convert-ignore' keeps whatever is left

            # 5. tee the sanitized record to disk (if requested)
            if out_fh is not None:
                out_fh.write(_format_record(header, seq))

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
        in_fh.close()
