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
import itertools
import os
import string
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

# Most characters of an offending line quoted in an error message. A file
# that is not FASTA at all (compressed, binary) can put megabytes on a
# single "line", and nobody wants that in a traceback.
_MAX_QUOTED_LINE_LENGTH = 60

# Upper-casing table for the 26 ASCII letters only. str.upper() applies the
# full Unicode case mapping, under which some non-ASCII characters turn into
# valid residues ('\u017f' -> 'S', '\u0131' -> 'I', '\u00df' -> 'SS', ...)
# and would then sail through validation. Sequences that contain anything
# beyond ASCII are upper-cased with this table instead, so those characters
# reach invalid_sequence_action exactly as they appear in the file.
_ASCII_UPPERCASE = str.maketrans(string.ascii_lowercase, string.ascii_uppercase)

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
    except (OSError, ValueError) as e:
        # ValueError is what open() raises for a path with a NUL byte in it
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
        When appending to an existing file whose last line has no line
        break, one is written first, so that the first appended header
        starts a line of its own rather than being glued onto the end of
        the file's last sequence (where it would be read back as sequence
        data, silently merging two records).  Default ``False``.

    Returns
    -------
    file object
        An open, writeable text file handle with a 1 MiB buffer.

    Raises
    ------
    ProtfastaException
        If the file cannot be opened (missing directory, permissions, ...).
    """
    needs_line_break = append and _ends_without_line_break(filename)
    try:
        fh = open(filename, 'a' if append else 'w', buffering=_WRITE_BUFFER, encoding=_WRITE_ENCODING, errors=_ENCODING_ERRORS)
    except (OSError, ValueError) as e:
        # ValueError is what open() raises for a path with a NUL byte in it
        raise ProtfastaException('Unable to open file for writing: %s\nException: %s' % (filename, e))

    if needs_line_break:
        fh.write('\n')
    return fh


####################################################################################################
#
#
def _ends_without_line_break(filename: PathLike) -> bool:
    """Return ``True`` if *filename* is a non-empty file whose last byte is not a line break.

    Only regular files are examined (reading a named pipe, say, would block),
    and any error while looking is treated as "no", leaving the subsequent
    open to report anything that is really wrong.

    Parameters
    ----------
    filename : str or os.PathLike
        The file about to be appended to.

    Returns
    -------
    bool
        Whether a line break must be written before appending.
    """
    try:
        if not os.path.isfile(filename):
            return False
        with open(filename, 'rb') as fh:
            if fh.seek(0, os.SEEK_END) == 0:
                return False
            fh.seek(-1, os.SEEK_END)
            return fh.read(1) not in (b'\n', b'\r')
    except OSError:
        return False


####################################################################################################
#
#
def _same_file(path_a: PathLike, path_b: PathLike) -> bool:
    """Return ``True`` if two paths refer to the same file on disk.

    Comparing the paths themselves - even fully resolved ones - is not
    enough: on a case-insensitive filesystem (the macOS and Windows
    defaults) ``Data.fasta`` and ``data.fasta`` are the same file, and a
    hard link gives one file two unrelated names.  This compares the files
    (device and inode) instead, and treats a path that does not exist as
    different from everything, since it cannot be an existing file.

    Parameters
    ----------
    path_a, path_b : str or os.PathLike
        The two paths.

    Returns
    -------
    bool
        Whether both paths name the same existing file.
    """
    try:
        return os.path.samefile(path_a, path_b)
    except (OSError, ValueError):
        # (ValueError: a path with a NUL byte in it, which cannot exist)
        return False


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
    empty_sequence_action: str = 'fail',
    leading_content_action: str = 'fail',
) -> None:
    """Validate the arguments shared by :func:`read_fasta` and :func:`read_fasta_stream`.

    This is a stateless guard function: it either returns ``None`` when
    every argument is acceptable, or raises a
    :class:`~protfasta.protfasta_exceptions.ProtfastaException` describing
    the first problem it encounters.

    If new functionality is added to ``read_fasta`` or
    ``read_fasta_stream``, the corresponding keyword must be validated
    here.

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

    empty_sequence_action : str, optional
        How to handle a header with no sequence after it.  Must be one
        of ``'ignore'``, ``'fail'``, or ``'remove'``.  ``'ignore'``
        cannot be combined with *output_filename*, because a record with
        an empty sequence cannot be written to a FASTA file.  Default
        ``'fail'``.

    leading_content_action : str, optional
        How to handle anything other than blank lines before the first
        header.  Must be ``'fail'`` or ``'remove'``.  Default ``'fail'``.

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

    # check the empty_sequence_action
    if empty_sequence_action not in ['ignore', 'fail', 'remove']:
        raise ProtfastaException("keyword 'empty_sequence_action' must be one of 'ignore','fail','remove'")

    # an empty sequence cannot be written out (write_fasta rejects it), so
    # catch this pairing up front rather than after the whole file has been
    # parsed
    if empty_sequence_action == 'ignore' and output_filename is not None:
        raise ProtfastaException("Cannot combine empty_sequence_action='ignore' with output_filename: records with no sequence cannot be written to a FASTA file")

    # check the leading_content_action. There is no 'ignore' here: content
    # with no header has nowhere to go in a list of [header, sequence] pairs,
    # so it can only be refused or thrown away
    if leading_content_action not in ['fail', 'remove']:
        raise ProtfastaException("keyword 'leading_content_action' must be one of 'fail','remove'")

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

        # alignment=True promises that gap characters come through untouched,
        # which a dictionary that converts them would quietly break
        if alignment is True:
            for key in correction_dictionary:
                if '-' in key:
                    raise ProtfastaException("Cannot combine alignment=True (which preserves '-' as a gap character) with a correction_dictionary that converts '-' (key %r)" % (key,))


####################################################################################################
#
#
def _upper_case_sequence(seq: str) -> str:
    """Upper-case the ASCII letters of a sequence, and nothing else.

    ``str.upper()`` is only safe on pure-ASCII text; see the note on
    ``_ASCII_UPPERCASE`` at the top of this module for what it does to
    the rest of Unicode.  ``str.isascii()`` is a constant-time flag
    check, so the overwhelmingly common pure-ASCII sequence still takes
    the fast ``upper()`` path.

    Parameters
    ----------
    seq : str
        A sequence as read from the file.

    Returns
    -------
    str
        The sequence with ``a``-``z`` replaced by ``A``-``Z``.  The
        length never changes.
    """
    if seq.isascii():
        return seq.upper()
    return seq.translate(_ASCII_UPPERCASE)


####################################################################################################
#
#
def _skip_leading_content(
    content: Iterable[str],
    leading_content_action: str = 'fail',
) -> tuple[Iterator[str], int]:
    """Deal with anything that comes before the first header in a FASTA file.

    A FASTA file should open with a ``">"`` header line.  If it doesn't,
    whatever comes first belongs to no record, and quietly throwing it
    away hides real problems: a file that is gzipped, UTF-16 or not FASTA
    at all has no line starting with ``">"``, so it would otherwise parse
    "successfully" into zero records, and a first header that has picked
    up a leading space would take its whole record with it.

    Shared by the load-it-all path (:func:`_parse_fasta_all`) and the
    streaming path (:func:`_stream_fasta`), which both run it immediately
    before :func:`_iter_fasta_lines` so the two always treat leading
    content the same way.  It uses the same definition of a header line
    as the parser: the first character is ``">"``.  Blank and
    whitespace-only lines are never counted as content.

    Parameters
    ----------
    content : Iterable[str]
        Iterable yielding the lines of a FASTA file.  It is consumed up
        to and including the first header line.

    leading_content_action : str, optional
        ``'fail'`` (default) raises on the first non-blank line that
        comes before a header.  ``'remove'`` discards every such line.

    Returns
    -------
    tuple[Iterator[str], int]
        The lines from the first header onwards (empty if the input has
        no header at all), and the number of non-blank lines discarded
        ahead of it.  The count is always ``0`` under ``'fail'``.

    Raises
    ------
    ProtfastaException
        If *leading_content_action* is ``'fail'`` and a non-blank line
        comes before the first header.
    """
    lines = iter(content)
    n_lines_discarded = 0

    for line_number, line in enumerate(lines, start=1):
        stripped = line.rstrip()

        if not stripped:
            continue

        if stripped[0] == '>':
            # Hand the header line back in front of everything not yet read.
            return itertools.chain((line,), lines), n_lines_discarded

        if leading_content_action == 'fail':
            # repr() rather than _printable(): this line may well be binary,
            # and repr() escapes control characters and lone surrogates and
            # makes leading whitespace visible.
            quoted = repr(stripped[:_MAX_QUOTED_LINE_LENGTH])
            if len(stripped) > _MAX_QUOTED_LINE_LENGTH:
                quoted = quoted + '...'
            raise ProtfastaException("Found content before the first header (line %i: %s). A FASTA file must start with a '>' header line, so this file may be compressed, in another encoding or format, or missing its first header. Pass leading_content_action='remove' to discard anything before the first header" % (line_number, quoted))

        n_lines_discarded += 1

    # No header anywhere, so there are no records to parse.
    return iter(()), n_lines_discarded


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
    lines.  Trailing whitespace is stripped from every line and blank
    lines are skipped.  Anything before the first header is ignored
    here; deciding whether such content is an error is the job of
    :func:`_skip_leading_content`, which both callers run first.

    A header that is not followed by any sequence data is yielded with
    an empty sequence (``''``).  The parser itself never drops a record;
    the callers decide what to do with empty ones via
    ``empty_sequence_action`` (see :func:`_keep_empty_record`).

    Parameters
    ----------
    content : Iterable[str]
        Iterable yielding the lines of a FASTA file.

    header_parser : callable or None, optional
        A function ``(str) -> str`` applied to every raw header (with the
        leading ``">"`` already removed).  ``None`` means headers are used
        verbatim.  The function must return a string for every header;
        anything else (``None`` from a regular expression that did not
        match, say) is an error rather than a silently dropped record.

    Yields
    ------
    tuple[str, str]
        ``(header, sequence)`` pairs in file order.  Sequences are
        concatenated from any multi-line runs and their ASCII letters
        upper-cased (see :func:`_upper_case_sequence`); a header
        with no sequence lines after it gives an empty sequence.

    Raises
    ------
    ProtfastaException
        If *header_parser* raises on, or returns something other than a
        string for, a header encountered in the file.
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
            # Flush the previous record. A header with no sequence lines
            # is yielded with an empty sequence rather than dropped, so
            # that the caller can report it (empty_sequence_action).
            if header is not None:
                yield (header, _upper_case_sequence(''.join(seq_parts)))

            # Start the new record.
            h = line[1:]
            if header_parser is not None:
                raw = h
                try:
                    h = header_parser(raw)
                except Exception as e:
                    raise ProtfastaException('header_parser raised an exception on header [%s]\nException: %s' % (_utilities._printable(raw), e))

                # A non-string return value (None from a failed regular
                # expression match is the classic case) must not be allowed
                # through: a None header would make the record look like
                # sequence data with no header and it would silently vanish.
                if not isinstance(h, str):
                    raise ProtfastaException('header_parser returned %s rather than a string for header [%s]' % (type(h).__name__, _utilities._printable(raw)))
            header = h
            seq_parts = []
        else:
            seq_parts.append(line)

    # Flush the final record (again, even if it has no sequence lines).
    if header is not None:
        yield (header, _upper_case_sequence(''.join(seq_parts)))


####################################################################################################
#
#
def _keep_empty_record(header: str, empty_sequence_action: str, emptied_by_conversion: bool = False) -> bool:
    """Decide what happens to a record with an empty sequence.

    Shared by the load-it-all path and the streaming path
    (:func:`_stream_fasta`) so the two always treat empty records the
    same way.  It is only called for records whose sequence is the empty
    string, which happens in two places: straight out of the parser, for
    a header with no sequence after it, and after the invalid-residue
    step, for a sequence made up entirely of characters the conversion
    table deletes (``'***'``, say).

    Parameters
    ----------
    header : str
        The header of the empty record (used in the error message).

    empty_sequence_action : str
        One of ``'fail'``, ``'remove'`` or ``'ignore'``.

    emptied_by_conversion : bool, optional
        ``True`` when the sequence was emptied by residue conversion
        rather than being empty in the file.  Only changes the wording
        of the error message.  Default ``False``.

    Returns
    -------
    bool
        ``True`` if the record should be kept (``'ignore'``), ``False``
        if it should be dropped (``'remove'``).

    Raises
    ------
    ProtfastaException
        If *empty_sequence_action* is ``'fail'``.
    """
    if empty_sequence_action == 'fail':
        if emptied_by_conversion:
            raise ProtfastaException("Converting invalid residues left no sequence at all for header (%s). Pass empty_sequence_action='remove' to skip such records, or 'ignore' to keep them with an empty sequence" % (_utilities._printable(header)))
        raise ProtfastaException("Found a header with no sequence after it (%s). Pass empty_sequence_action='remove' to skip such records, or 'ignore' to keep them with an empty sequence" % (_utilities._printable(header)))

    return empty_sequence_action == 'ignore'


####################################################################################################
#
#
def _parse_fasta_all(
    content: Iterable[str],
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    verbose: bool = False,
    empty_sequence_action: str = 'fail',
    leading_content_action: str = 'fail',
) -> list[list[str]]:
    """Parse FASTA content into a list of ``[header, sequence]`` pairs.

    A thin wrapper over :func:`_iter_fasta_lines` that deals with anything
    before the first header (:func:`_skip_leading_content`), materialises
    every record and, when asked, checks header uniqueness as it goes.  Used by
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

    empty_sequence_action : str, optional
        What to do with a header that has no sequence after it:
        ``'fail'`` (default) raises, ``'remove'`` drops the record and
        ``'ignore'`` keeps it with an empty sequence.  Empty records are
        dealt with before the header-uniqueness check, so a removed
        record never counts towards a duplicate header.

    leading_content_action : str, optional
        What to do with non-blank lines before the first header:
        ``'fail'`` (default) raises and ``'remove'`` discards them.

    Returns
    -------
    list[list[str]]
        A list of two-element lists ``[header, sequence]``.

    Raises
    ------
    ProtfastaException
        If *expect_unique_header* is ``True`` and a duplicate header is
        found, if *empty_sequence_action* is ``'fail'`` and a header
        has no sequence, or if *leading_content_action* is ``'fail'``
        and there is content before the first header.
    """

    return_data: list[list[str]] = []
    n_empty_removed = 0

    # anything before the first header (raises under 'fail')
    content, n_leading_lines_removed = _skip_leading_content(content, leading_content_action)

    # Only allocate a header-tracking set when we actually need it; the set
    # holds references to the header strings already in return_data, so it
    # costs one pointer-sized slot per record rather than a copy.
    seen_headers: Optional[set[str]] = set() if expect_unique_header else None

    # Every [header, sequence] pair is a GC-tracked container, so a large
    # file makes the cyclic garbage collector run thousands of times to scan
    # objects that cannot possibly be part of a cycle - about 20% of the
    # parse time at ten million records. The loop below is bounded and
    # creates no cycles of its own, so the collector is paused for its
    # duration and restored (to whatever state it was in) on the way out,
    # exceptions included. (A header_parser runs inside the loop, but any
    # cyclic garbage it did create would simply wait for the next
    # collection after the parse rather than being lost.) The streaming
    # parser does not do this, because arbitrary amounts of user code run
    # between its yields.
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        for header, seq in _iter_fasta_lines(content, header_parser):
            # a header with no sequence after it (raises under 'fail')
            if not seq and not _keep_empty_record(header, empty_sequence_action):
                n_empty_removed += 1
                continue

            if seen_headers is not None:
                if header in seen_headers:
                    raise ProtfastaException('Found duplicate header (%s)' % (_utilities._printable(header)))
                seen_headers.add(header)
            return_data.append([header, seq])
    finally:
        if gc_was_enabled:
            gc.enable()

    if verbose:
        if leading_content_action == 'remove':
            print('[INFO]: Removed %i lines of content before the first header' % (n_leading_lines_removed))
        if empty_sequence_action == 'remove':
            print('[INFO]: Removed %i records with no sequence' % (n_empty_removed))
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
    empty_sequence_action: str = 'fail',
    leading_content_action: str = 'fail',
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

    empty_sequence_action : str, optional
        What to do with a header that has no sequence after it; see
        :func:`_parse_fasta_all`.  Default ``'fail'``.

    leading_content_action : str, optional
        What to do with non-blank lines before the first header; see
        :func:`_parse_fasta_all`.  Default ``'fail'``.

    Returns
    -------
    list[list[str]]
        A list of two-element lists ``[header, sequence]`` in the order
        they appear in the file.  Sequences are upper-cased.

    Raises
    ------
    ProtfastaException
        If the file cannot be opened, a duplicate header is detected
        (when *expect_unique_header* is ``True``), a header has no
        sequence (when *empty_sequence_action* is ``'fail'``), or there
        is content before the first header (when
        *leading_content_action* is ``'fail'``).
    """

    fh = _open_fasta(filename)

    if verbose:
        print('[INFO]: Read in file %s (streaming)' % (filename))

    with fh:
        return _parse_fasta_all(fh,
                                expect_unique_header=expect_unique_header,
                                header_parser=header_parser,
                                verbose=verbose,
                                empty_sequence_action=empty_sequence_action,
                                leading_content_action=leading_content_action)


####################################################################################################
#
#
def _iter_fasta(
    filename: PathLike,
    header_parser: Optional[Callable[[str], str]] = None,
) -> Iterator[tuple[str, str]]:
    """Yield raw ``(header, sequence)`` pairs from a FASTA file, streaming.

    A convenience wrapper that opens *filename* and hands it to
    :func:`_iter_fasta_lines`.  No duplicate detection, empty-record
    handling, invalid-residue handling, or alignment-gap logic is
    performed -- each record is yielded exactly as parsed (so a header
    with no sequence gives an empty sequence).  The public streaming entry point,
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
def _check_encodable(text: str, description: str) -> None:
    """Raise if *text* cannot be written to a FASTA file as UTF-8.

    Output is encoded with the ``surrogateescape`` handler, so the lone
    surrogates that stand in for non-UTF-8 bytes read from a file go back
    out as those bytes.  Any *other* lone surrogate (``'\\ud800'``, say)
    has no encoding at all, and would otherwise blow up part-way through
    the write with the output file already truncated.

    Parameters
    ----------
    text : str
        A header or a sequence.

    description : str
        How to refer to *text* in the error message.

    Raises
    ------
    ProtfastaException
        If *text* cannot be encoded.
    """
    try:
        text.encode(_WRITE_ENCODING, _ENCODING_ERRORS)
    except UnicodeEncodeError:
        raise ProtfastaException('%s contains a character that cannot be written as UTF-8' % (description))


####################################################################################################
#
#
def _check_wrapped_lines_survive(header: str, seq: str, linelength: int) -> None:
    """Raise if wrapping *seq* would change what the file means when read back.

    Two things can go wrong once a sequence is cut into lines: a line
    that starts with ``">"`` is read back as a header, splitting the
    record in two, and whitespace at the end of a line is stripped by
    the parser, shortening the sequence.  Neither can happen to a
    validated sequence, but ``invalid_sequence_action='ignore'`` and a
    direct :func:`protfasta.write_fasta` call both let anything through.

    Only the first and last character of each output line matter, and
    those sit at fixed strides through the sequence, so they are pulled
    out with two slices rather than by building every line.

    Parameters
    ----------
    header : str
        The record header (used in the error message).

    seq : str
        The sequence about to be written.  Must not be empty.

    linelength : int
        Number of residues per output line; any falsy value means the
        sequence is written on a single line.

    Raises
    ------
    ProtfastaException
        If any output line would start with ``">"`` or end in whitespace.
    """
    # a '>' anywhere at all is rare, so one scan settles it for nearly
    # every sequence without looking at where the lines start
    if '>' in seq:
        line_starts = seq[::linelength] if linelength else seq[0]
        if '>' in line_starts:
            raise ProtfastaException("Sequence associated with [%s] cannot be written: one of its lines would start with '>' and be read back as a header" % (_utilities._printable(header)))

    # the last character of every full line, plus the last character of the
    # sequence (which ends the final, possibly shorter, line)
    if linelength:
        line_ends = seq[linelength - 1::linelength] + seq[-1]
    else:
        line_ends = seq[-1]

    # Letters are not whitespace, so isalpha() clears the usual case in one
    # short pass. Otherwise (a gap or a '*' at a line end, say) look properly:
    # split() with no argument splits on any whitespace, so it only hands
    # the string back whole when there is none in it.
    if not line_ends.isalpha():
        if line_ends.split() != [line_ends]:
            raise ProtfastaException('Sequence associated with [%s] cannot be written: one of its lines would end in whitespace, which is stripped when the file is read back' % (_utilities._printable(header)))


####################################################################################################
#
#
def _check_record_is_writable(header: object, seq: object, linelength: int = 60) -> None:
    """Raise unless a record can be written and read back unchanged.

    Shared by :func:`protfasta.write_fasta` and the streaming tee in
    :func:`_stream_fasta`, so the two can never disagree about what is
    fit to go into a FASTA file.  The header and sequence must be
    strings, the sequence must contain something other than whitespace
    (the reader skips blank lines, so a blank sequence would leave a
    header with no sequence), neither may contain a line break (which
    would be read back as a record boundary), both must be encodable,
    and no line of the wrapped sequence may start with ``">"`` or end in
    whitespace (see :func:`_check_wrapped_lines_survive`).

    A clean record costs a handful of C-level operations: ``isspace()``
    stops at the first non-whitespace character, ``isascii()`` is a flag
    lookup (so nothing is test-encoded unless it has to be), and the
    wrapping check looks only at the first and last character of each
    output line.

    Parameters
    ----------
    header : object
        The record header.  Anything but a string is rejected.

    seq : object
        The sequence.  Anything but a string is rejected.

    linelength : int, optional
        Number of residues per output line, as passed to
        :func:`_format_record`.  Default ``60``.

    Raises
    ------
    ProtfastaException
        If the record fails any of the checks above.
    """
    printable = _utilities._printable

    if not isinstance(header, str):
        raise ProtfastaException('Header [%s] is not a string (got %s)' % (header, type(header).__name__))
    if '\n' in header or '\r' in header:
        raise ProtfastaException('Header [%s] contains a line break' % (printable(header)))
    if not isinstance(seq, str):
        raise ProtfastaException('Sequence associated with [%s] is not a string (got %s)' % (printable(header), type(seq).__name__))
    if not seq or seq.isspace():
        raise ProtfastaException('Sequence associated with [%s] is empty (or only whitespace)' % (printable(header)))
    if '\n' in seq or '\r' in seq:
        raise ProtfastaException('Sequence associated with [%s] contains a line break' % (printable(header)))

    if not header.isascii():
        _check_encodable(header, 'Header [%s]' % (printable(header)))
    if not seq.isascii():
        _check_encodable(seq, 'Sequence associated with [%s]' % (printable(header)))

    _check_wrapped_lines_survive(header, seq, linelength)


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
    empty_sequence_action: str = 'fail',
    leading_content_action: str = 'fail',
) -> Iterator[Union[tuple[str, str], list[str]]]:
    """Stream a FASTA file record-by-record with full sanitization.

    This is the streaming engine behind :func:`protfasta.read_fasta_stream`.
    It applies the same processing steps as :func:`protfasta.read_fasta`
    -- header-uniqueness checks, duplicate-record and duplicate-sequence
    handling, and invalid-residue handling -- but yields one record at a
    time instead of materializing the whole dataset.

    The processing order matches :func:`protfasta.read_fasta`:

    0. Content before the first header (*leading_content_action*).
    1. Records with no sequence (*empty_sequence_action*).
    2. Header uniqueness (*expect_unique_header*).
    3. Duplicate records (*duplicate_record_action*).
    4. Duplicate sequences (*duplicate_sequence_action*).
    5. Invalid residues (*invalid_sequence_action*).
    6. Optional tee to *output_filename*.

    Peak memory is ``O(number of records)`` for the auxiliary
    duplicate/uniqueness bookkeeping (16-byte digests -- never whole
    sequences -- plus the header strings themselves under
    *expect_unique_header* and ``duplicate_sequence_action='fail'``, which
    both need to name a header later) and ``O(single record)`` for the
    sequence data itself.  When *expect_unique_header* is ``False`` and
    all duplicate actions are ``'ignore'``, the auxiliary bookkeeping is
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

    empty_sequence_action : str, optional
        What to do with a header that has no sequence after it:
        ``'fail'`` (default) raises at that record, ``'remove'`` drops it
        and ``'ignore'`` yields it with an empty sequence.  This is a
        per-record decision, so it needs no extra memory.

    leading_content_action : str, optional
        What to do with non-blank lines before the first header:
        ``'fail'`` (default) raises before any record is yielded and
        ``'remove'`` discards them.  Needs no extra memory either.

    Yields
    ------
    tuple[str, str] or list[str]
        ``(header, sequence)`` pairs (or ``[header, sequence]`` lists when
        *return_list* is ``True``) in file order.  Sequences are
        upper-cased and sanitized.

    Raises
    ------
    ProtfastaException
        On content before the first header (for
        ``leading_content_action='fail'``),
        a header with no sequence (for ``empty_sequence_action='fail'``),
        a duplicate header/record/sequence (for the relevant ``'fail'``
        actions) or an invalid residue (for ``'fail'``/``'convert'``).
        Because parsing is lazy, these are raised mid-iteration, at the
        offending record.
    """

    # Only allocate bookkeeping structures for the actions that need them.
    seen_headers: Optional[set[str]] = set() if expect_unique_header else None

    # duplicate records: a set of 16-byte digests of (header, sequence).
    # (Unlike read_fasta, the stream cannot compare records directly, since
    # it does not keep the records it has already yielded.) When headers are
    # required to be unique a duplicate record (same header AND sequence) is
    # impossible - the header check fires first - so the record check would
    # only ever hash every sequence for nothing.
    seen_records: Optional[set[bytes]] = (
        set() if (duplicate_record_action in ('fail', 'remove') and not expect_unique_header) else None
    )

    # duplicate sequences. For 'fail' the first header seen for each digest
    # is retained so the error message can name both offenders; for 'remove'
    # only the digests are needed, and not holding on to the headers leaves
    # every header string free to be released once it has been yielded.
    seq_lookup: Optional[dict[bytes, str]] = {} if duplicate_sequence_action == 'fail' else None
    seen_seqs: Optional[set[bytes]] = set() if duplicate_sequence_action == 'remove' else None

    # Build the converter once (it caches its translate table) rather than
    # once per record.
    convert = _utilities._make_converter(correction_dictionary, alignment) if invalid_sequence_action.startswith('convert') else None
    check_valid = _utilities.check_sequence_is_valid
    printable = _utilities._printable

    n_read = 0
    n_yielded = 0
    n_empty_removed = 0
    n_dup_records_removed = 0
    n_dup_seqs_removed = 0
    n_invalid_removed = 0
    n_converted = 0

    # Open the input first so a missing input file never creates an empty
    # output file.
    in_fh = _open_fasta(filename)
    out_fh: Optional[IO[str]] = None

    try:
        # 0. anything before the first header (raises under 'fail'). Done
        # before the output is opened, so an input that turns out not to be
        # a FASTA file never creates or truncates the output file.
        lines, n_leading_lines_removed = _skip_leading_content(in_fh, leading_content_action)

        if output_filename is not None:
            out_fh = _open_output(output_filename)

        if verbose:
            print('[INFO]: Streaming file %s' % (filename))

        for header, seq in _iter_fasta_lines(lines, header_parser):
            n_read += 1

            # 1. a header with no sequence after it (raises under 'fail').
            # Dealt with first, as in read_fasta, so a removed empty record
            # never counts towards a duplicate header.
            if not seq and not _keep_empty_record(header, empty_sequence_action):
                n_empty_removed += 1
                continue

            # 2. header uniqueness
            if seen_headers is not None:
                if header in seen_headers:
                    raise ProtfastaException('Found duplicate header (%s)' % (printable(header)))
                seen_headers.add(header)

            # 3. duplicate records (identical header AND sequence)
            if seen_records is not None:
                key = _utilities._record_hash(header, seq)
                if key in seen_records:
                    if duplicate_record_action == 'fail':
                        raise ProtfastaException('Found duplicate entries of the following record:\n>%s\n%s' % (printable(header), printable(seq)))
                    n_dup_records_removed += 1
                    continue
                seen_records.add(key)

            # 4. duplicate sequences (identical sequence, any header)
            if seq_lookup is not None:
                digest = _utilities._seq_hash(seq)
                if digest in seq_lookup:
                    raise ProtfastaException('Found duplicate sequences associated with the following headers\n1. %s\n\n2. %s' % (printable(seq_lookup[digest]), printable(header)))
                seq_lookup[digest] = header
            elif seen_seqs is not None:
                digest = _utilities._seq_hash(seq)
                if digest in seen_seqs:
                    n_dup_seqs_removed += 1
                    continue
                seen_seqs.add(digest)

            # 5. invalid-residue handling (per record)
            if invalid_sequence_action == 'ignore':
                pass

            elif invalid_sequence_action == 'fail':
                (status, info) = check_valid(seq, alignment)
                if status is not True:
                    raise ProtfastaException('Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (_utilities._describe_character(str(info)), printable(header), printable(seq)))

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
                        inner = 'Failed on invalid amino acid: %s\nTaken from entry...\n>%s\n%s\n' % (_utilities._describe_character(str(info)), printable(header), printable(seq))
                        raise ProtfastaException("\n\n******* Despite fixing fixable errors, additional problems remain with the sequence*********\n%s" % (inner))

                elif invalid_sequence_action == 'convert-remove':
                    (status, _info) = check_valid(seq, alignment)
                    if status is not True:
                        n_invalid_removed += 1
                        continue

                # 'convert-ignore' keeps whatever is left

                # conversion deletes characters, so it can leave nothing
                # behind (raises under 'fail')
                if not seq and not _keep_empty_record(header, empty_sequence_action, emptied_by_conversion=True):
                    n_empty_removed += 1
                    continue

            # 6. tee the sanitized record to disk (if requested), under the
            # same checks write_fasta applies, so the stream can never write
            # a file that write_fasta would have refused
            if out_fh is not None:
                _check_record_is_writable(header, seq)
                out_fh.write(_format_record(header, seq))

            n_yielded += 1
            if return_list:
                yield [header, seq]
            else:
                yield (header, seq)

        if verbose:
            if leading_content_action == 'remove':
                print('[INFO]: Removed %i lines of content before the first header' % (n_leading_lines_removed))
            if empty_sequence_action == 'remove':
                print('[INFO]: Removed %i of %i records with no sequence' % (n_empty_removed, n_read))
            if duplicate_record_action == 'remove':
                print('[INFO]: Removed %i of %i due to duplicate records ' % (n_dup_records_removed, n_read))
            if duplicate_sequence_action == 'remove':
                print('[INFO]: Removed %i of %i due to duplicate sequences ' % (n_dup_seqs_removed, n_read))
            if invalid_sequence_action in ('convert', 'convert-ignore', 'convert-remove'):
                print('[INFO]: Converted invalid residues in %i sequences' % (n_converted))
            if invalid_sequence_action in ('remove', 'convert-remove'):
                print('[INFO]: Removed %i of %i due to sequences with invalid characters' % (n_invalid_removed, n_read))
            print('[INFO]: Streamed %i of %i records from %s' % (n_yielded, n_read, filename))

    finally:
        if out_fh is not None:
            out_fh.close()
        in_fh.close()
