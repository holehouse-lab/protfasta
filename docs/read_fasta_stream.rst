read_fasta_stream
=================

``read_fasta_stream`` is the **streaming** counterpart to
:func:`protfasta.read_fasta`. It takes the same set of arguments, but
instead of reading the whole file and returning a fully-materialised
``dict`` or ``list``, it returns a **generator** that yields one
sanitized record at a time. This keeps peak memory bounded and makes it
possible to process FASTA files far larger than RAM (tens or hundreds of
millions of sequences).

At its simplest::

    import protfasta

    for header, sequence in protfasta.read_fasta_stream('huge.fasta'):
        # process each record on the fly
        if len(sequence) > 1000:
            print(header, len(sequence))

Because a lazily-produced result cannot be a dictionary, the
``dict``-versus-``list`` distinction of ``read_fasta`` collapses:
``read_fasta_stream`` yields ``(header, sequence)`` tuples by default, or
``[header, sequence]`` lists when ``return_list=True``.


read_fasta vs read_fasta_stream
................................

The two functions take the same arguments (``read_fasta_stream`` adds
one more, ``silence_warnings``, and uses different defaults for the
duplicate checks - see below) and apply the same sanitization steps in
the same order. Choose between them based on the size of the file and
the shape of result you need:

    *  :func:`protfasta.read_fasta` loads the file and returns every
       record at once as a ``dict`` or ``list``. Use it whenever the
       dataset fits comfortably in memory - it is the more convenient
       interface and lets you index and re-iterate the result freely.
    *  ``read_fasta_stream`` yields records one at a time and never
       holds more than a single record's sequence data in memory. Use
       it for files that do not fit in memory, or when you only need a
       single pass over the records.


What differs under streaming
.............................

The full parameter set of :func:`protfasta.read_fasta` is supported, but
a few arguments behave slightly differently because a stream never sees
the whole file up front:

    *  **Duplicate-checking defaults.** ``read_fasta_stream`` defaults to
       ``expect_unique_header=False`` and
       ``duplicate_record_action='ignore'``, whereas ``read_fasta``
       defaults to ``expect_unique_header=True`` and
       ``duplicate_record_action='fail'``. **Duplicate headers and
       duplicate records are therefore not detected by default when
       streaming.** This is deliberate - each of those checks costs
       ``O(records)`` memory, which would defeat the point of a streamer
       built for files larger than RAM. See `Memory characteristics`_
       below.
    *  **Return type.** ``return_list=False`` (default) yields
       ``(header, sequence)`` tuples; ``return_list=True`` yields
       ``[header, sequence]`` lists. There is no dictionary return type.
    *  **When errors are raised.** All *argument* validation is eager -
       an invalid keyword combination, a ``filename`` that does not
       exist, or an ``output_filename`` that is the input file, raises
       immediately, before any record is produced. *Data-dependent*
       failures, however, are inherent to streaming and are raised
       **mid-iteration**, at the offending record: a header with no
       sequence under ``empty_sequence_action='fail'``, a duplicate header, a
       duplicate record or sequence under a ``'fail'`` action, an invalid
       residue under ``invalid_sequence_action='fail'`` / ``'convert'``,
       or a ``header_parser`` that raises on (or does not return a string
       for) a header. Records yielded before the failure have already
       been handed to you, and written to ``output_filename`` if one was
       given.
    *  **verbose output.** Per-total summaries (for example, "removed 5
       of 100 duplicate records") can only be reported once the
       generator has been fully consumed, so they are emitted at
       exhaustion rather than up front.
    *  **output_filename.** When provided, each sanitized record is
       teed to disk as it is yielded. The output file is therefore only
       complete once the generator has been fully consumed, and it must
       be a different file from the input ``filename`` (opening it would
       otherwise truncate the input before it had been read). The check
       compares the files themselves rather than their paths, so a
       symlink, a hard link or, on a case-insensitive filesystem such as
       the macOS default, a differently-capitalised name for the input
       is rejected too. The input is opened before the
       output is created, so a problem opening the input never leaves an
       empty output file behind. The output is byte-for-byte what
       :func:`protfasta.write_fasta` would produce with its default
       ``linelength``.


Memory characteristics
......................

With the default arguments, ``read_fasta_stream`` is **flat in memory**.
Only one record is held at a time and no per-record bookkeeping is kept,
so peak memory is independent of file size and files far larger than RAM
stream without issue.

The duplicate- and uniqueness-handling checks are still available, but
each one has to remember what it has already seen, so enabling any of
them adds auxiliary state that grows with the *number of records*:

    *  ``expect_unique_header=True`` keeps a running set of every header;
    *  ``duplicate_record_action`` set to ``'fail'`` or ``'remove'``
       keeps a running set of 16-byte record digests. (When
       ``expect_unique_header=True`` this check is skipped entirely,
       since unique headers already rule out duplicate records.)
    *  ``duplicate_sequence_action`` set to ``'fail'`` or ``'remove'``
       keeps a running set of 16-byte sequence digests (``'fail'`` also
       keeps the first header seen for each sequence, so that its error
       can name both records).

No sequence is ever stored, only digests and (where noted) headers, so
this is still far lighter than a full load: the benchmarks in the
repository's ``benchmark/`` directory measure roughly 120 bytes per record
for ``duplicate_sequence_action='remove'`` and 200 bytes per record for
``expect_unique_header=True`` with UniProt-style headers, against about 800
bytes per record for :func:`protfasta.read_fasta`. It does grow with the
file, though - 1-2 GB per ten million records.
Whenever one of these checks is enabled a one-time warning is emitted
naming the responsible keyword(s); pass ``silence_warnings=True`` to
suppress it. The check itself is always performed either way.

Note that ``expect_unique_header=True`` on its own is sufficient - it
transparently promotes the default ``duplicate_record_action='ignore'``
to ``'fail'``, since unique headers already preclude duplicate records.

To restore full ``read_fasta``-equivalent checking while streaming:

.. code-block:: python

    import protfasta

    stream = protfasta.read_fasta_stream('input.fasta',
                                         expect_unique_header=True)


Basic usage
...............

.. code-block:: python

    import protfasta

    for header, sequence in protfasta.read_fasta_stream('huge.fasta'):
        if len(sequence) > 1000:
            print(header, len(sequence))


Using a header parser
......................

As with ``read_fasta``, the optional ``header_parser`` argument is a
callable ``(str) -> str`` applied to every raw header. This is commonly
used to extract a UniProt accession from a structured header:

.. code-block:: python

    import protfasta

    def uniprot_id(header):
        # '>sp|P12345|NAME_HUMAN ...' -> 'P12345'
        return header.split('|')[1] if '|' in header else header

    for acc, seq in protfasta.read_fasta_stream('uniprot.fasta',
                                                header_parser=uniprot_id):
        ...

The fallback matters: the parser is smoke-tested with a plain string
before the file is opened (see ``check_header_parser``), so a version
that assumes a ``|`` is present would fail that test. The parser must
also return a string for every header - returning ``None`` (from a
regular expression that did not match, say) raises a
``ProtfastaException`` naming the header rather than silently dropping
the record.


Sanitizing on the fly
......................

Every sanitization option available to ``read_fasta`` works while
streaming. For example, to stream a very large file while converting
non-standard residues and dropping duplicate sequences:

.. code-block:: python

    import protfasta

    stream = protfasta.read_fasta_stream('huge.fasta',
                                         invalid_sequence_action='convert',
                                         duplicate_sequence_action='remove')

    for header, seq in stream:
        ...


Streaming filter + write
.........................

You can stream through an input file, keep only the records you care
about, and write them out, holding only the records you keep in memory
rather than the whole file:

.. code-block:: python

    import protfasta

    keep = []
    for header, seq in protfasta.read_fasta_stream('huge.fasta'):
        if 100 <= len(seq) <= 500:
            keep.append([header, seq])

    protfasta.write_fasta(keep, 'filtered.fasta')

Alternatively, let ``read_fasta_stream`` write the sanitized output for
you as it goes via ``output_filename`` - just remember the file is only
complete once the generator has been fully consumed:

.. code-block:: python

    import protfasta

    stream = protfasta.read_fasta_stream('huge.fasta',
                                         invalid_sequence_action='convert',
                                         output_filename='clean.fasta')

    # fully consume the stream so the output file is written in full
    for header, seq in stream:
        ...


For usage examples see the :doc:`examples` page. Full API documentation
is shown below.


Documentation
...............

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodule:: protfasta
   :noindex:

.. autofunction:: read_fasta_stream
