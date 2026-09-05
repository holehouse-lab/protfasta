from __future__ import annotations

import os
import warnings
from typing import Callable, Iterable, Iterator, Optional, Union

from protfasta import utilities as _utilities
from protfasta import io as _io
from protfasta._configs import STANDARD_AAS, STANDARD_CONVERSION
from protfasta import protfasta as _protfasta
from protfasta.protfasta_exceptions import ProtfastaException

_ROOT = os.path.abspath(os.path.dirname(__file__))


## ------------------------------------------------------------
# READTHEDOCS versioning hack
#
# Read the Docs builds the documentation from a bare checkout without
# installing the package, so the versioningit-generated _version.py does not
# exist there. In that environment (and only there) the version is computed
# from git on the fly instead. Both paths are anchored on this file's
# location rather than the working directory so the check does not depend on
# where sphinx happens to be run from.
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile(os.path.join(_ROOT, '_version.py')):
    import versioningit
    __version__ = versioningit.get_version(os.path.dirname(_ROOT))
else:
    from ._version import __version__


## ------------------------------------------------------------


__all__ = [
    'read_fasta',
    'read_fasta_stream',
    'write_fasta',
    'ProtfastaException',
    'STANDARD_AAS',
    'STANDARD_CONVERSION',
    '__version__',
]


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
    filename: _io.PathLike,
    expect_unique_header: bool = True,
    header_parser: Optional[Callable[[str], str]] = None,
    check_header_parser: bool = True,
    duplicate_sequence_action: str = 'ignore',
    duplicate_record_action: str = 'fail',
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    return_list: bool = False,
    output_filename: Optional[_io.PathLike] = None,
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

    Files are decoded as UTF-8 regardless of the platform locale.  A
    leading byte-order mark is stripped, and bytes that are not valid
    UTF-8 are carried through (rather than raising a decode error) so
    that a stray byte in a sequence is reported as an invalid residue
    like any other and a non-UTF-8 header survives a read/write round
    trip byte-for-byte.

    Parameters
    ----------
    filename : str or os.PathLike
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
        begins.  The function must return a string for every header in
        the file; if it raises, or returns anything else (``None`` from
        a regular expression that did not match, say), a
        ``ProtfastaException`` naming the offending header is raised.

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

    output_filename : str, os.PathLike, or None, optional
        If provided, the final (sanitized) set of sequences is written
        to a new FASTA file at this path before the function returns.

    correction_dictionary : dict or None, optional
        A mapping of non-standard characters to replacement strings used
        when *invalid_sequence_action* involves conversion.  When
        ``None``, the built-in table is used:

        * ``B`` -> ``N``, ``U`` -> ``C``, ``X`` -> ``G``, ``Z`` -> ``Q``
        * ``*`` -> ``''``, ``-`` -> ``''``, ``' '`` -> ``''``

        A custom dictionary **replaces** the built-in table entirely
        (an empty dictionary is treated the same as ``None``).  Keys
        must be non-empty strings and values must be strings.

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
    _io.check_filename(filename)
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

    # first deal with duplicate records. If every header is unique (which
    # internal_parse_fasta_file has just enforced) a duplicate record - same
    # header AND same sequence - cannot exist, so the pass would only hash
    # every sequence to find nothing; skip it.
    if expect_unique_header:
        updated = raw
        if verbose and duplicate_record_action == 'remove':
            print('[INFO]: Removed 0 of %i due to duplicate records (headers are unique)' % (len(raw)))
    else:
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
        

    # If we wanted to write the final set of sequences we're going to use...
    # (tested against None rather than for truthiness so that an empty
    # string is reported as an unwritable path instead of silently skipped)
    if output_filename is not None:
        write_fasta(updated, output_filename)

    # if we asked for a list...
    if return_list is True:
        return updated

    return _utilities.convert_list_to_dictionary(updated, verbose)



# ------------------------------------------------------------------
#
def read_fasta_stream(
    filename: _io.PathLike,
    expect_unique_header: bool = False,
    header_parser: Optional[Callable[[str], str]] = None,
    check_header_parser: bool = True,
    duplicate_sequence_action: str = 'ignore',
    duplicate_record_action: str = 'ignore',
    invalid_sequence_action: str = 'fail',
    alignment: bool = False,
    return_list: bool = False,
    output_filename: Optional[_io.PathLike] = None,
    correction_dictionary: Optional[dict[str, str]] = None,
    verbose: bool = False,
    silence_warnings: bool = False,
) -> Iterator[Union[tuple[str, str], list[str]]]:
    """Stream a FASTA file record-by-record, sanitizing as it goes.

    This is the streaming counterpart to :func:`read_fasta`.  It takes the
    same arguments but, instead of loading the whole file and returning a
    dictionary or list, it returns a **generator** that yields one sanitized
    record at a time.  This keeps peak memory bounded and makes it possible to
    process files far larger than RAM::

        for header, seq in read_fasta_stream('huge.fasta'):
            ...

    .. warning::

       **The duplicate/uniqueness-checking defaults differ from**
       :func:`read_fasta`.  To keep streaming flat in memory by default,
       this function defaults to ``expect_unique_header=False`` and
       ``duplicate_record_action='ignore'``, whereas :func:`read_fasta`
       defaults to ``expect_unique_header=True`` and
       ``duplicate_record_action='fail'``.  As a result, **duplicate headers
       and duplicate records are not detected by default when streaming.**
       This is deliberate: those checks each cost ``O(records)`` memory (a
       running set of headers / digests), which would defeat the point of a
       streamer built for files larger than RAM.  Re-enable them explicitly
       (e.g. ``expect_unique_header=True``) if you need them -- a one-time
       warning will then flag the memory cost.  See the Memory section below.

    Because a lazily-produced result cannot be a dictionary, the
    dict-versus-list distinction of :func:`read_fasta` collapses: this
    function yields ``(header, sequence)`` tuples by default, or
    ``[header, sequence]`` lists when *return_list* is ``True``.

    Memory
    ------
    **By default this function is flat in memory**: only one record is held
    at a time and no per-record bookkeeping is kept, so peak memory is
    independent of file size and files far larger than RAM stream fine.  This
    is why the duplicate/uniqueness-checking defaults differ from
    :func:`read_fasta` (see the warning above).

    The *duplicate/uniqueness checks* are still available, but each needs to
    remember what it has already seen, so enabling one adds ``O(records)``
    auxiliary state (hundreds of MB or more on very large files):

    * ``expect_unique_header=True`` keeps a running set of every header;
    * ``duplicate_record_action`` in ``('fail', 'remove')`` keeps a
      running set of 16-byte record digests (skipped entirely when
      ``expect_unique_header=True``, since unique headers already rule
      out duplicate records);
    * ``duplicate_sequence_action`` in ``('fail', 'remove')`` keeps a
      running set of sequence digests.

    Whenever a memory-growing check is enabled a one-time warning is emitted
    (suppress it with ``silence_warnings=True``).

    All argument validation is performed eagerly, at call time, so bad
    keyword combinations (and a *filename* that does not exist) raise
    immediately, before iteration begins.  Data-dependent errors -- a
    duplicate header, a duplicate record or sequence under a ``'fail'``
    action, or an invalid residue -- are inherent to streaming and are
    therefore raised **mid-iteration**, at the offending record, rather
    than up front.  The input file is opened before *output_filename* is
    created, so a failure to open the input never leaves an empty output
    file behind.

    The same processing steps as :func:`read_fasta` are applied, in the
    same order (header uniqueness, duplicate records, duplicate sequences,
    invalid residues).  See :func:`read_fasta` for a full description of
    each argument; the notes below cover only where streaming differs.

    Parameters
    ----------
    filename : str or os.PathLike
        Path to the FASTA file to read.

    expect_unique_header : bool, optional
        As in :func:`read_fasta`, but **defaults to ``False`` here** so that
        streaming is flat in memory by default.  When ``True`` a running set
        of seen headers is kept (``O(records)`` memory), and the default
        ``duplicate_record_action='ignore'`` is promoted to ``'fail'`` (unique
        headers already preclude duplicate records), so ``True`` works on its
        own.  Default ``False``.

    header_parser : callable or None, optional
        As in :func:`read_fasta`.  Default ``None``.

    check_header_parser : bool, optional
        As in :func:`read_fasta`.  Default ``True``.

    duplicate_sequence_action : str, optional
        As in :func:`read_fasta`.  The ``'fail'`` and ``'remove'``
        variants keep a running set of 16-byte sequence digests (never
        whole sequences), i.e. ``O(records)`` memory.  Default ``'ignore'``.

    duplicate_record_action : str, optional
        As in :func:`read_fasta`, but **defaults to ``'ignore'`` here** so
        that streaming is flat in memory by default.  The ``'fail'`` and
        ``'remove'`` variants keep a running set of 16-byte record
        digests (``O(records)`` memory).  Default ``'ignore'``.

    invalid_sequence_action : str, optional
        As in :func:`read_fasta`.  Every mode is a per-record decision
        and streams without extra state.  Default ``'fail'``.

    alignment : bool, optional
        As in :func:`read_fasta`.  Default ``False``.

    return_list : bool, optional
        If ``True``, yield ``[header, sequence]`` lists (matching the
        internal format of :func:`read_fasta` with ``return_list=True``);
        otherwise yield ``(header, sequence)`` tuples.  Default
        ``False``.

    output_filename : str, os.PathLike, or None, optional
        If provided, each sanitized record is written to this path as it
        is yielded (60 residues per line, as in :func:`write_fasta`).
        The file is only complete once the generator has been fully
        consumed.  Must differ from *filename*.

    correction_dictionary : dict or None, optional
        As in :func:`read_fasta`.  Default ``None``.

    verbose : bool, optional
        If ``True``, print an opening message and, once the generator is
        exhausted, a summary of removed/converted counts.  Per-total
        summaries are only emitted at exhaustion because a stream never
        sees the file total up front.  Default ``False``.

    silence_warnings : bool, optional
        If ``True``, suppress the one-time warning emitted when a
        memory-growing duplicate/uniqueness check is enabled (see the
        Memory section above).  Default ``False``.

    Returns
    -------
    Iterator[tuple[str, str]] or Iterator[list[str]]
        A generator yielding ``(header, sequence)`` tuples (or
        ``[header, sequence]`` lists when *return_list* is ``True``) in
        file order.  Sequences are upper-cased and sanitized.

    Raises
    ------
    ProtfastaException
        Eagerly, for invalid argument combinations; lazily (during
        iteration) for data-dependent failures.

    See Also
    --------
    read_fasta : Load and return the whole file as a dict or list.
    write_fasta : Write sequences out to a FASTA file.
    """

    # Now that the streaming defaults are flat (duplicate_record_action='ignore'),
    # a caller who only flips expect_unique_header=True would otherwise trip the
    # 'cannot expect unique headers and ignore duplicate records' guard below.
    # Since unique headers already preclude duplicate records, transparently
    # promote the default 'ignore' -> 'fail' so expect_unique_header=True works on
    # its own.
    if expect_unique_header and duplicate_record_action == 'ignore':
        duplicate_record_action = 'fail'

    # Validate arguments eagerly (this function is not itself a generator,
    # so its body runs at call time -- bad keywords fail fast, before any
    # iteration begins, matching read_fasta).
    _io.check_filename(filename)
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

    # A missing input is an argument problem, not a data problem, so report
    # it now rather than at the first next() call. Anything else that can go
    # wrong at open time (a directory, permissions) is still caught, but
    # lazily, by _open_fasta.
    if not os.path.exists(filename):
        raise ProtfastaException('Unable to find file: %s' % (filename))

    # Streaming-specific guard: reading and simultaneously overwriting the
    # same file would corrupt the input mid-stream.  read_fasta is immune
    # because it reads the whole file before writing, but streaming is not.
    # realpath (rather than abspath) so a symlink to the input is caught too.
    if output_filename is not None:
        if os.path.realpath(output_filename) == os.path.realpath(filename):
            raise ProtfastaException("keyword 'output_filename' must differ from 'filename' when streaming")

    # Warn (once, eagerly) if a memory-growing check is enabled, so the caller
    # knows their peak memory will scale with the number of records rather than
    # staying flat. The checks are still performed -- this only surfaces the cost
    # and the flat-memory recipe.
    if not silence_warnings:
        growing = []
        if expect_unique_header:
            growing.append("expect_unique_header=True")
        elif duplicate_record_action in ('fail', 'remove'):
            # (with unique headers a duplicate record is impossible, so the
            # record check is skipped and costs nothing extra)
            growing.append("duplicate_record_action=%r" % duplicate_record_action)
        if duplicate_sequence_action in ('fail', 'remove'):
            growing.append("duplicate_sequence_action=%r" % duplicate_sequence_action)
        if growing:
            warnings.warn(
                "read_fasta_stream(): %s require O(number of records) bookkeeping, "
                "so peak memory grows with file size rather than staying flat. For "
                "truly flat, file-size-independent streaming set expect_unique_header="
                "False, duplicate_record_action='ignore' and duplicate_sequence_action="
                "'ignore'. Pass silence_warnings=True to suppress this message."
                % ", ".join(growing),
                stacklevel=2,
            )

    return _io._stream_fasta(filename,
                             expect_unique_header=expect_unique_header,
                             header_parser=header_parser,
                             duplicate_sequence_action=duplicate_sequence_action,
                             duplicate_record_action=duplicate_record_action,
                             invalid_sequence_action=invalid_sequence_action,
                             alignment=alignment,
                             return_list=return_list,
                             output_filename=output_filename,
                             correction_dictionary=correction_dictionary,
                             verbose=verbose)



# ------------------------------------------------------------------
#
def write_fasta(
    fasta_data: Union[dict[str, str], list[list[str]]],
    filename: _io.PathLike,
    linelength: Union[int, bool, None] = 60,
    append_to_fasta: bool = False,
) -> None:
    """Write sequences to a FASTA file.

    Accepts sequence data as either a dictionary (header -> sequence) or
    a list of ``[header, sequence]`` pairs and writes a standards-
    compliant FASTA file to *filename*.

    The data is validated in full **before** the file is opened, so a bad
    entry raises without creating, truncating, or appending to
    *filename*.  Every header and sequence must be a string, every
    sequence must be non-empty, and neither may contain a line break (a
    break inside either would be read back as a record boundary, silently
    corrupting the file).  Output is written as UTF-8; headers that came
    from a non-UTF-8 file via :func:`read_fasta` are written back out
    with their original bytes.

    Parameters
    ----------
    fasta_data : dict[str, str] or list[list[str]]
        Sequence data.  If a dictionary, keys are headers and values are
        amino-acid sequences.  If a list, each element must be a
        two-element list (or tuple) ``[header, sequence]``.

    filename : str or os.PathLike
        Destination file path.  Should conventionally end with
        ``.fasta`` or ``.fa``, but this is not enforced.

    linelength : int, bool, or None, optional
        Maximum number of residues per line in the output.  Default is
        ``60`` (the UniProt convention).  Values from ``1`` to ``4`` are
        clamped to ``5``.  Set to ``0`` (or any negative value),
        ``None``, or ``False`` to write each sequence on a single line.
        A numerical string such as ``'60'`` is accepted and cast.

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
        If *fasta_data* is neither a dictionary nor a list, if a list
        element is not a two-item list or tuple, if a header or sequence
        is not a string or contains a line break, if a sequence is
        empty, if *linelength* is not an integer (or
        ``0``/``None``/``False``), if *append_to_fasta* is not a
        boolean, or if *filename* cannot be opened for writing.
    """

    # This part of code means we can pass either a dictionary or a list of lists
    # in for write_fasta to deal with. Either way, `records` iterates over
    # (header, sequence) pairs.
    records: Iterable[tuple[str, str]]
    if isinstance(fasta_data, dict):
        records = fasta_data.items()

    elif isinstance(fasta_data, list):
        # quick validate: every element must be a [header, sequence] pair (a
        # tuple is fine too). The type matters as well as the length - a bare
        # two-character string would otherwise unpack into a one-character
        # header and a one-character sequence and be written out as such.
        for i in fasta_data:
            if not isinstance(i, (list, tuple)) or len(i) != 2:
                raise ProtfastaException('While processing a list for write_fasta_file at least one of the elements was not a 2-position sublist:\n%s' % (str(i)))
        records = fasta_data  # type: ignore[assignment]

    else:
        raise ProtfastaException("keyword 'fasta_data' must be a dictionary of header:sequence pairs or a list of [header, sequence] pairs (got %s)" % (type(fasta_data).__name__))

    if not isinstance(append_to_fasta, bool):
        raise ProtfastaException("keyword 'append_to_fasta' must be a boolean")

    # override line length for sane input. Note the int() cast happens before
    # any numerical comparison so that a non-numerical linelength raises a
    # protfasta exception rather than an opaque TypeError.
    if linelength is None or linelength is False:
        linelength = False

    else:
        # cast linelength to an integer here as a soft type checking, and
        # if it's shorter than 5 then reset to 5
        try:
            linelength = int(linelength)
        except (TypeError, ValueError):
            raise ProtfastaException("keyword 'linelength' must be an integer, or one of 0/None/False (got %s)" % (repr(linelength)))

        if linelength < 1:
            linelength = False
        elif linelength < 5:
            linelength = 5

    # Validate every entry before touching the filesystem, so a bad entry
    # never leaves a truncated file (or a half-appended one) behind. Headers
    # and sequences must be strings (anything else either crashes inside the
    # formatting below or, worse, gets str()-ed into the file), sequences
    # must be non-empty, and neither may contain a line break, which would be
    # read back as a record boundary.
    printable = _utilities._printable
    for header, seq in records:
        if not isinstance(header, str):
            raise ProtfastaException('Header [%s] is not a string (got %s)' % (header, type(header).__name__))
        if '\n' in header or '\r' in header:
            raise ProtfastaException('Header [%s] contains a line break' % (printable(header)))
        if not isinstance(seq, str):
            raise ProtfastaException('Sequence associated with [%s] is not a string (got %s)' % (printable(header), type(seq).__name__))
        if not seq:
            raise ProtfastaException('Sequence associated with [%s] is empty' % (printable(header)))
        if '\n' in seq or '\r' in seq:
            raise ProtfastaException('Sequence associated with [%s] contains a line break' % (printable(header)))

    # If append_to_fasta is False the file is truncated/created; otherwise
    # new entries are added to the end of an existing file.
    with _io._open_output(filename, append=append_to_fasta) as fh:
        write = fh.write
        for header, seq in records:
            # Each record is assembled into a single string (header line,
            # line-length-sliced sequence, blank separator) and emitted with
            # one write call - a handful of C-level operations per record
            # rather than a Python-level write per line.
            write(_io._format_record(header, seq, linelength))
    
