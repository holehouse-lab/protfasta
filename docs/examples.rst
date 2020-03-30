.. _example_label:

Examples
=================


.. toctree::
   :maxdepth: 2
   :caption: Contents:


It's often easies to see how to use code through some well-worked examples. Here we provide some simple examples that illustrate how **protfasta** can be used to read and write FASTA files.


``read_fasta`` examples
..........

Some possible exampls of reading FASTA files using protfasta:



**Example 1: Simple read in FASTA file**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta')


**Example 2: Simple read in FASTA file and ignore duplicate FASTA records and return a nested-list of residues**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta', 
                                     expect_unique_header=False, 
				     return_list=True, 
				     duplicate_record_action='ignore')


**Example 3: Read in FASTA file and correct invalid residues using standard error correction dictionary**

.. code-block:: python

    import protfasta

    sequences = protfasta.read_fasta('inputfile.fasta', 
                                     invalid_sequence_action='convert')


**Example 4: Read in FASTA file and correct invalid residues using a custom dictionary**

.. code-block:: python

    import protfasta

    CD = {'U':'G', '-':''}
    sequences = protfasta.read_fasta('inputfile.fasta', 
                                     invalid_sequence_action='convert', 
				     correction_dictionary=CD)



``write_fasta`` examples
.............


.. code-block:: python

    # input example using a sequence dictionary
    import protfasta

    sequence_in = {'seq1': 'MEEPQSDPSVEPPLS', 'seq2': 'DEAPRMPEAAPPVAPA'}
    protfasta.write_fasta(sequence_in, 'example.fasta')

.. code-block:: python    
	   
    # input example using a sequence list
    import protfasta

    sequence_in = [['seq1','MEEPQSDPSVEPPLS'], ['seq2', 'DEAPRMPEAAPPVAPA']]
    protfasta.write_fasta(sequence_in, 'example.fasta')
