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


``iter_fasta`` examples (large files)
.......................................

For files that are too large to hold in memory, use
``iter_fasta`` to stream records one at a time.

**Example 10 - Stream a huge FASTA file**

.. code-block:: python

    import protfasta

    long_seq_count = 0
    for header, seq in protfasta.iter_fasta('metagenome.fasta'):
        if len(seq) > 1000:
            long_seq_count += 1
    print(long_seq_count)


**Example 11 - Streaming filter + write**

.. code-block:: python

    import protfasta

    out = []
    for header, seq in protfasta.iter_fasta('huge.fasta'):
        if 100 <= len(seq) <= 500:
            out.append([header, seq])

    protfasta.write_fasta(out, 'filtered.fasta')


``write_fasta`` examples
...........................

**Example 12 - Write from a dictionary**

.. code-block:: python

    import protfasta

    sequence_in = {'seq1': 'MEEPQSDPSVEPPLS',
                   'seq2': 'DEAPRMPEAAPPVAPA'}
    protfasta.write_fasta(sequence_in, 'example.fasta')


**Example 13 - Write from a list**

.. code-block:: python

    import protfasta

    sequence_in = [['seq1', 'MEEPQSDPSVEPPLS'],
                   ['seq2', 'DEAPRMPEAAPPVAPA']]
    protfasta.write_fasta(sequence_in, 'example.fasta')


**Example 14 - Single-line sequences**

.. code-block:: python

    import protfasta

    protfasta.write_fasta(sequence_in, 'example.fasta', linelength=None)


**Example 15 - Append to an existing FASTA file**

.. code-block:: python

    import protfasta

    protfasta.write_fasta(sequence_in, 'archive.fasta',
                          append_to_fasta=True)
