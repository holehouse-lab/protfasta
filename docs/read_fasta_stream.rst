read_fasta_stream
=================

``read_fasta_stream`` is the **streaming** counterpart to
:func:`protfasta.read_fasta`. It takes an identical set of arguments,
but instead of reading the whole file and returning a fully-materialised
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

The two functions share the same signature and apply the same
sanitization steps in the same order. Choose between them based on the
size of the file and the shape of result you need:

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

    *  **Return type.** ``return_list=False`` (default) yields
       ``(header, sequence)`` tuples; ``return_list=True`` yields
       ``[header, sequence]`` lists. There is no dictionary return type.
    *  **When errors are raised.** All *argument* validation is eager -
       an invalid keyword combination raises immediately, before any
       record is produced. *Data-dependent* failures, however, are
       inherent to streaming and are raised **mid-iteration**, at the
       offending record: a duplicate header, a duplicate record or
       sequence under a ``'fail'`` action, or an invalid residue under
       ``invalid_sequence_action='fail'`` / ``'convert'``.
    *  **verbose output.** Per-total summaries (for example, "removed 5
       of 100 duplicate records") can only be reported once the
       generator has been fully consumed, so they are emitted at
       exhaustion rather than up front.
    *  **output_filename.** When provided, each sanitized record is
       teed to disk as it is yielded. The output file is therefore only
       complete once the generator has been fully consumed, and it must
       differ from the input ``filename``.


Memory characteristics
......................

Peak memory is ``O(single record)`` for the sequence data itself. The
duplicate- and uniqueness-handling actions add auxiliary bookkeeping
that grows with the *number of records* (headers plus 16-byte sequence
digests - never whole sequences), so they are still far lighter than a
full load. When ``expect_unique_header=False`` and every duplicate
action is ``'ignore'``, no bookkeeping is allocated at all and memory is
flat regardless of file size.


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
        return header.split('|')[1]

    for acc, seq in protfasta.read_fasta_stream('uniprot.fasta',
                                                header_parser=uniprot_id):
        ...


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
about, and write them out in bounded memory:

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
