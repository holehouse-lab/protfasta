.. _example_label:

Examples
=================


.. toctree::
   :maxdepth: 2
   :caption: Contents:


It is often easier to see how to use code through some well-worked
examples. Here we provide worked examples for every major feature of
**protfasta**.


``read_fasta`` examples
...........................

**Example 1 - Simple read**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta')
    # sequences is a dict: {header: sequence}


**Example 2 - Allow duplicate FASTA records and return a list**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta',
                                     expect_unique_header=False,
                                     return_list=True,
                                     duplicate_record_action='ignore')


**Example 3 - Correct invalid residues using the standard table**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta',
                                     invalid_sequence_action='convert')


**Example 4 - Correct invalid residues using a custom dictionary**

.. code-block:: python

    import protfasta

    CD = {'U': 'G', '-': ''}
    sequences = protfasta.read_fasta('inputfile.fasta',
                                     invalid_sequence_action='convert',
                                     correction_dictionary=CD)


**Example 5 - Convert what can be converted, drop the rest**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta',
                                     invalid_sequence_action='convert-remove')


**Example 6 - Parse an alignment (keep gap characters)**

.. code-block:: python

    import protfasta

    aln = protfasta.read_fasta('alignment.fasta',
                               alignment=True,
                               invalid_sequence_action='convert')


**Example 7 - Custom header parser**

Extract just the UniProt accession from a structured header such as
``>sp|P12345|NAME_HUMAN ...``:

.. code-block:: python

    import protfasta

    def get_accession(header):
        return header.split('|')[1]

    sequences = protfasta.read_fasta('uniprot.fasta',
                                     header_parser=get_accession)


**Example 8 - Remove duplicate sequences and write directly to disk**

.. code-block:: python

    import protfasta

    protfasta.read_fasta('inputfile.fasta',
                         duplicate_sequence_action='remove',
                         output_filename='unique.fasta')


**Example 9 - Fast read with no sanity checking**

By default **protfasta** performs a lot of sanity checking. When you
already trust the input file, you can turn it all off:

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta',
                                     invalid_sequence_action='ignore',
                                     duplicate_record_action='ignore',
                                     duplicate_sequence_action='ignore',
                                     expect_unique_header=False)


``read_fasta_stream`` examples (large files)
.............................................

``read_fasta_stream`` is the streaming counterpart to ``read_fasta``. It
takes the same arguments and applies the same sanitization, but instead
of returning a ``dict`` / ``list`` it returns a generator that yields one
``(header, sequence)`` record at a time (or ``[header, sequence]`` lists
with ``return_list=True``).

**When should I use** ``read_fasta_stream`` **instead of** ``read_fasta``\ **?**

    *  Reach for ``read_fasta`` by default. It loads the whole file and
       hands you every record at once, so you can index it, re-iterate
       it, and pass it around freely. This is the right choice whenever
       the dataset fits comfortably in memory.
    *  Reach for ``read_fasta_stream`` when the file is too large to hold
       in memory, or when a **single pass** over the records is all you
       need (counting, filtering, transforming, or re-writing). Peak
       memory stays bounded to roughly one record regardless of file
       size.

A useful rule of thumb: if you would immediately loop over the result of
``read_fasta`` and never keep the whole collection, ``read_fasta_stream``
does the same job without ever materialising it.

**Example 10 - Stream a huge FASTA file**

.. code-block:: python

    import protfasta

    # Count long sequences in a file far larger than RAM without ever
    # holding more than one record in memory at a time.
    long_seq_count = 0
    for header, seq in protfasta.read_fasta_stream('metagenome.fasta'):
        if len(seq) > 1000:
            long_seq_count += 1
    print(long_seq_count)


**Example 11 - Streaming filter + write**

.. code-block:: python

    import protfasta

    out = []
    for header, seq in protfasta.read_fasta_stream('huge.fasta'):
        if 100 <= len(seq) <= 500:
            out.append([header, seq])

    protfasta.write_fasta(out, 'filtered.fasta')


**Example 12 - Sanitize while streaming**

Every sanitization option available to ``read_fasta`` works while
streaming. Here we convert non-standard residues and drop duplicate
sequences as records flow past:

.. code-block:: python

    import protfasta

    stream = protfasta.read_fasta_stream('huge.fasta',
                                         invalid_sequence_action='convert',
                                         duplicate_sequence_action='remove')

    for header, seq in stream:
        # 'seq' is already converted; duplicates have been skipped
        ...


**Example 13 - Let read_fasta_stream write the cleaned file for you**

Passing ``output_filename`` tees each sanitized record to disk as it is
yielded. Because the write happens lazily, you must fully consume the
generator for the output file to be complete:

.. code-block:: python

    import protfasta

    stream = protfasta.read_fasta_stream('huge.fasta',
                                         invalid_sequence_action='convert',
                                         output_filename='clean.fasta')

    # Exhaust the stream so every record is written out.
    for header, seq in stream:
        ...


``write_fasta`` examples
...........................

**Example 14 - Write from a dictionary**

.. code-block:: python

    import protfasta

    sequence_in = {'seq1': 'MEEPQSDPSVEPPLS',
                   'seq2': 'DEAPRMPEAAPPVAPA'}
    protfasta.write_fasta(sequence_in, 'example.fasta')


**Example 15 - Write from a list**

.. code-block:: python

    import protfasta

    sequence_in = [['seq1', 'MEEPQSDPSVEPPLS'],
                   ['seq2', 'DEAPRMPEAAPPVAPA']]
    protfasta.write_fasta(sequence_in, 'example.fasta')


**Example 16 - Single-line sequences**

.. code-block:: python

    import protfasta

    protfasta.write_fasta(sequence_in, 'example.fasta', linelength=None)


**Example 17 - Append to an existing FASTA file**

.. code-block:: python

    import protfasta

    protfasta.write_fasta(sequence_in, 'archive.fasta',
                          append_to_fasta=True)
