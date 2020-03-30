pfasta
=================

**pfasta** is a command-line tool for working with FASTA files to filter and sanitize them based on various criterion. This includes:

    *  Filtering out sequences that contain invalid amino acids
    *  Take sequences that contain invalid characters and replace/fix them
    *  Filter a set of sequences by a maximum and/or minimum sequence length
    *  Sub-sample a set of sequences for building a reduced set of randomly selected sequences

At it's basline, **pfasta** takes a single sequence file as input and writes a new output sequence. There are a series of flags that can be applied, as outlined in the Usage section below. 

Usage
.........

.. code-block:: python

    pfasta <flags> filename.fasta




.. code-block:: none

    -o <output filename> (default: output.fasta)
       Define the name of the output FASAT file

    --non-unique-header
      Flag that, if provided allows multiple FASTA records to have identical headers

    --duplicate-record (default: fail)
      Flag that provides a keyword that defines how duplicate FASTA records are dealt with. 
      Options are:
          fail   : throws an exception and exits the parsing 
	  ignore : duplicate records are retained
	  remove : duplicate records are removed

    --duplicate-sequence (default: fail)
      Flag that provides a keyword that defines how duplicate sequences are dealt with. 
      Options are:
          fail   : throws an exception and exits the parsing 
	  ignore : duplicate sequences are retained
	  remove : duplicate sequences are removed

    --invalid-sequence (default: fail)
      Flag that provides a keyword that defines how invalid sequences are dealt with. 
      Options are:
          fail                : throws an exception and exits the parsing 
	  ignore              : invalid sequences are retained
	  remove              : invalid sequences are removed
	  convert-all         : invalid residues are converted according to the standard conversion table 
	                        (shown below) but if OTHER invalid residues are found an exception is raised
	                        B->N,    U->C,    X->G,    Z->Q,    '*'->'',    '-'->''
	  convert-res         : invalid residues are converted according to the standard conversion table
	                        with the exception of sequence-alignment gaps ('-') 
	  convert-all-ignore  : invalid residues are converted according to the standard conversion table,
                                and if OTHER invalid residues are found they are ignored
	  convert-res-ignore  : invalid residues are converted according to the standard conversion table,
	                        with the exception of the sequence-aligment gap ('-') character, but 
			        if OTHER invalid residues are found they are ignored

    --number-lines (default: 60)
      Flag that defines the number of lines in the output FASTA file

    --shortest-seq-lines (default: None)
      Flag that defines a filter that sets the shortest sequence returned

    --longest-seq-lines (default: None)
      Flag that defines a filter that sets the longest sequence returned

    --random-subsample (default: None)
      Flag that defines the number of randomly sub-sampled sequences. Allows a test FASTA file to be 
      generated as a sub-set for testing analysis pipelines

    --print-statistics
      Flag that, if provided, means statistics about the FINAL set of sequences written

    --no-outputfile
      Flag that, if provided, means NO outputfile is generated.

    --silent
      Flag that, if provided, means pfasta generates ZERO output to STDOUT 




.. toctree::
   :maxdepth: 2
   :caption: Contents:

		     
