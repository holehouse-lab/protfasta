read_fasta
=================

``read_fasta`` is a one-stop-shop for reading in FASTA files! Customizable keywords allow a variety of sanitizing functions which include:

    *  Ignore, remove, or convert sequences with invalid amino acid characters (``B``/``U``/``X``/``*``/``-``)
    *  Ignore or remove duplicate sequences or duplicate FASTA records
    *  Alternatively, allow duplicate sequences, headers, and FASTA records (something most other parsers do not)
    *  Arbitrary conversion of amino acids via a customizable ``correction_dictionary``

Once parsed, ``read_fasta`` returns either a dictionary of header-to-sequence values or a nested list, where each sub-list contains two elements (header, sequence).

For usage examples see the :ref:`example_label` page. Full documentation is shown below. 


Documentation
...............

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodule:: protfasta

.. autofunction:: read_fasta
		     
