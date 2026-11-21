read_fasta
=================

``read_fasta`` is the primary entry point to **protfasta**. It reads a
FASTA file, optionally sanitizes its contents, and returns the records
either as a dictionary (``header -> sequence``) or as a list of
``[header, sequence]`` pairs.

At its simplest::

    import protfasta
    sequences = protfasta.read_fasta('proteins.fasta')

Many optional keyword arguments customize the behaviour - including
duplicate handling, invalid-residue correction, alignment-gap support,
custom header parsing, and automatic writing of the sanitized output.


What ``read_fasta`` can do
...........................

    *  Ignore, remove, convert, or fail on sequences containing
       non-standard amino-acid characters
       (``B``, ``U``, ``X``, ``Z``, ``*``, ``-``, ``' '``, ...).
    *  Ignore, remove, or fail on duplicate FASTA records
       (same header **and** same sequence).
    *  Ignore, remove, or fail on duplicate sequences (same sequence,
       different headers).
    *  Preserve alignment gap characters (``-``) when
       ``alignment=True``.
    *  Apply a caller-supplied ``header_parser`` function to every
       raw header (useful for extracting accession IDs).
    *  Optionally write the sanitized result to a new FASTA file via
       ``output_filename``.
    *  Override the built-in invalid-character conversion table via
       a custom ``correction_dictionary``.


Processing pipeline
.....................

Sanitization happens in a fixed order:

    1. File is streamed from disk, anything before the first header is
       processed according to ``leading_content_action``, headers are
       parsed with ``header_parser`` (if provided), records with no
       sequence are processed according to ``empty_sequence_action``,
       and header uniqueness is checked (when
       ``expect_unique_header=True``).
    2. Duplicate **records** are processed according to
       ``duplicate_record_action``.
    3. Duplicate **sequences** are processed according to
       ``duplicate_sequence_action``.
    4. **Invalid residues** are processed according to
       ``invalid_sequence_action``. If conversion leaves a sequence
       with nothing in it (``***``, say), that record is then handled
       by ``empty_sequence_action`` like any other empty one.
    5. The sanitized set is optionally written to
       ``output_filename``.
    6. The result is returned as a ``dict`` (default) or a ``list``
       of ``[header, sequence]`` pairs (when ``return_list=True``).

Two consequences of that order are worth knowing:

    *  The duplicate checks run on sequences as they are in the file,
       *before* any conversion. Two sequences that only become identical
       once converted (``ACB`` and ``ACN``, with ``B`` -> ``N``) are not
       duplicates as far as steps 2 and 3 are concerned, and both come
       back. If that matters, run the output through ``read_fasta``
       again.
    *  With ``expect_unique_header=False`` and the default dictionary
       return, records that share a header overwrite one another (the
       last one wins) and ``read_fasta`` emits a warning saying how many
       were lost. Pass ``return_list=True`` to keep them all.

Incompatible option combinations (for example,
``expect_unique_header=True`` together with
``duplicate_record_action='ignore'``) are caught before the file is
read.

Both ``filename`` and ``output_filename`` accept either a string or a
:class:`pathlib.Path`. Anything that goes wrong - a bad keyword, a
missing or unreadable file, an output file that cannot be created, an
invalid residue under ``invalid_sequence_action='fail'``, a header with
no sequence under ``empty_sequence_action='fail'``, content before the
first header under ``leading_content_action='fail'``, a
``header_parser`` that raises on (or returns something other than a
string for) a real header - is raised as a ``ProtfastaException``, so
callers only need to catch one exception type.


How the file is parsed
.......................

The same rules apply to :func:`protfasta.read_fasta` and
:func:`protfasta.read_fasta_stream`, which share a single parser:

    *  A line whose first character is ``>`` starts a new record; the
       rest of that line (after ``header_parser``, if given) is the
       header. A ``>`` anywhere else - including after leading
       whitespace - is treated as sequence data.
    *  All following lines up to the next header are the sequence. They
       are concatenated, so sequences may be wrapped at any width, and
       the letters ``a``-``z`` are upper-cased. Nothing else is touched:
       a non-ASCII character is never case-mapped into a valid residue
       (Python's ``upper()`` would turn ``ß`` into ``SS``), so it shows
       up as the invalid residue it is.
    *  Trailing whitespace is stripped from every line, and blank or
       whitespace-only lines are skipped. Whitespace *inside* a sequence
       line is kept, and is then an invalid residue like any other
       (spaces are removed by the conversion table below; tabs are not).
    *  Unix (``\n``), Windows (``\r\n``) and old Mac (``\r``) line
       endings are all accepted.
    *  Blank lines before the first header are fine. Anything else up
       there is handled by ``leading_content_action``: ``'fail'`` (the
       default) raises a ``ProtfastaException`` quoting the offending
       line, and ``'remove'`` throws it away. This is what catches a
       file that isn't really FASTA - gzipped, UTF-16, or just missing
       its first header - which has no line starting with ``>`` and so
       would otherwise read as zero records without a word. It also
       catches a first header that's picked up a leading space. Before
       version 0.1.27 this content was always dropped silently; pass
       ``leading_content_action='remove'`` to get that behaviour back.
    *  A header with no sequence lines after it (for example ``>a``
       followed directly by ``>b``, or a header on the last line of the
       file) is handled by ``empty_sequence_action``: ``'fail'`` (the
       default) raises a ``ProtfastaException`` naming the header,
       ``'remove'`` drops the record, and ``'ignore'`` keeps it with an
       empty sequence (this cannot be combined with ``output_filename``,
       since an empty sequence cannot be written to a FASTA file).
       Before version 0.1.25 such records were always dropped silently;
       pass ``empty_sequence_action='remove'`` to get that behaviour
       back.


Header parsing
...............

``header_parser`` is a callable ``(str) -> str`` applied to every raw
header (with the leading ``>`` already removed) before any uniqueness
check. It is smoke-tested with a plain string before the file is opened
(disable that with ``check_header_parser=False``), so a parser that
assumes a particular structure needs a fallback:

.. code-block:: python

    def get_accession(header):
        # '>sp|P12345|NAME_HUMAN ...' -> 'P12345'
        return header.split('|')[1] if '|' in header else header

The parser must return a string for every header in the file. Returning
anything else - ``None`` from a regular expression that did not match is
the classic case - raises a ``ProtfastaException`` naming the header
rather than silently dropping the record.


Encoding
.........

Files are always decoded as UTF-8, regardless of the platform locale
(so behaviour is identical on Linux, macOS and Windows), and written
back out as UTF-8. Two details make this robust to real-world files:

    *  A leading byte-order mark (which some Windows editors add) is
       stripped. Without this the BOM would hide the first ``>`` and the
       first record would silently vanish.
    *  Bytes that are not valid UTF-8 (a Latin-1 accented character in a
       header, say) do not raise a decode error. They are carried through
       and written back out as the original bytes, so a read/write round
       trip is lossless. A stray byte inside a *sequence* simply shows up
       as an invalid residue and is handled by ``invalid_sequence_action``
       like any other.


Default conversion table
..........................

When ``invalid_sequence_action`` includes conversion and no custom
``correction_dictionary`` is supplied, these replacements are applied
(a custom dictionary replaces this table entirely, although an empty one
is treated the same as ``None``; its keys must be non-empty strings and
its values strings, and multi-character keys are allowed):

    *  ``B`` -> ``N``
    *  ``U`` -> ``C``
    *  ``X`` -> ``G``
    *  ``Z`` -> ``Q``
    *  ``*`` -> ``''`` (removed)
    *  ``-`` -> ``''`` (removed; preserved if ``alignment=True``)
    *  ``' '`` -> ``''`` (spaces removed; other whitespace, such as a
       tab, is not in the table)

Two things to know about a custom ``correction_dictionary``:

    *  All the replacements happen in a single pass, so the output of
       one is never picked up by another: ``{'B': 'X', 'X': 'G'}`` turns
       ``B`` into ``X``, not ``G``, and the order of the dictionary
       doesn't matter. Where two keys match at the same spot (``X`` and
       ``XX``), the longer one wins. Before version 0.1.27 a dictionary
       with a multi-character key was applied one key after another
       instead, so its replacements could chain.
    *  With ``alignment=True`` it can't contain a key with ``-`` in it.
       Preserving gaps and converting them are contradictory, so that
       combination raises rather than picking one.


Large files
............

For files that do not fit comfortably in memory, consider using the
streaming parser :func:`protfasta.read_fasta_stream` instead.
``read_fasta`` itself streams the file from disk (so it will not load
the entire file as a single string), but it still builds an in-memory
data structure of all records; ``read_fasta_stream`` avoids that by
yielding one sanitized record at a time.


For usage examples see the :doc:`examples` page. Full API
documentation is shown below.


Documentation
...............

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodule:: protfasta

.. autofunction:: read_fasta
