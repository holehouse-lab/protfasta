iter_fasta
=================

``iter_fasta`` is a memory-efficient **streaming** parser for very large
FASTA files (tens or hundreds of millions of sequences) where loading
the entire dataset into memory with :func:`protfasta.read_fasta` is not
feasible.

Unlike ``read_fasta``, ``iter_fasta``:

    *  Yields ``(header, sequence)`` tuples one at a time instead of
       returning a fully materialised ``dict`` / ``list``.
    *  Keeps peak memory usage to O(single record), regardless of how
       large the input file is.
    *  Performs **no** duplicate detection, invalid-residue handling,
       or alignment-gap logic. Callers who need those features should
       apply them record-by-record as they consume the iterator.

If you are working with a FASTA file that fits comfortably in memory,
prefer :func:`protfasta.read_fasta` - it offers full sanitization and
is only marginally slower per record.


Basic usage
...............

.. code-block:: python

    import protfasta

    for header, sequence in protfasta.iter_fasta('huge.fasta'):
        # process each sequence on the fly
        if len(sequence) > 1000:
            print(header, len(sequence))


Using a header parser
......................

The optional ``header_parser`` argument is a callable ``(str) -> str``
applied to every raw header before it is yielded. This is commonly
used to extract a UniProt accession from a structured header:

.. code-block:: python

    import protfasta

    def uniprot_id(header):
        # '>sp|P12345|NAME_HUMAN ...' -> 'P12345'
        return header.split('|')[1]

    for acc, seq in protfasta.iter_fasta('uniprot.fasta',
                                         header_parser=uniprot_id):
        ...


Streaming filtering + writing
...............................

Because ``write_fasta`` accepts a list of ``[header, sequence]`` pairs,
you can stream through an input file, keep only the records you care
about, and write them out in bounded memory:

.. code-block:: python

    import protfasta

    keep = []
    for header, seq in protfasta.iter_fasta('huge.fasta'):
        if 100 <= len(seq) <= 500:
            keep.append([header, seq])
            # optionally flush every N records to bound memory further

    protfasta.write_fasta(keep, 'filtered.fasta')


Documentation
...............

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodule:: protfasta
   :noindex:

.. autofunction:: iter_fasta
