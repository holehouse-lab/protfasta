pfasta
=================

**pfasta** is a command-line tool for working with FASTA files to filter and
sanitize them based on various criteria. It is installed automatically with
**protfasta** and can be invoked as ``pfasta`` from the command line.

At its simplest, **pfasta** takes a single FASTA file as input and writes a
sanitized FASTA file as output. It can:

    *  Filter out (or convert, or fail on) sequences containing non-standard
       amino-acid characters
    *  Remove, ignore, or fail on duplicate FASTA records and duplicate
       sequences
    *  Filter sequences by minimum and/or maximum length
    *  Randomly sub-sample a set of sequences (useful for building a small
       test set from a large FASTA file)
    *  Print summary statistics (count, median / quartile / min / max length)
    *  Replace commas in FASTA headers with semicolons (helpful when the
       downstream pipeline treats the header as part of a CSV)

Usage
.........

.. code-block:: none

    pfasta <flags> filename.fasta


Command-line options
.....................

.. code-block:: none

    filename
        Positional argument: path to the input FASTA file.

    -o <output filename>                    (default: output.fasta)
        Output FASTA file.

    --non-unique-header
        If set, multiple FASTA records are allowed to share the same header.
        By default duplicate headers cause pfasta to fail.

    --duplicate-record {ignore,fail,remove} (default: fail)
        How to deal with duplicate records (same header AND same sequence):
            fail   - raise an exception and exit
            ignore - keep all duplicate records
            remove - keep only the first occurrence

    --duplicate-sequence {ignore,fail,remove} (default: ignore)
        How to deal with duplicate sequences (same sequence, any header):
            fail   - raise an exception and exit
            ignore - keep all duplicate sequences
            remove - keep only the first occurrence of each sequence

    --invalid-sequence <mode>                (default: fail)
        How to deal with non-standard amino-acid characters. Available
        modes:

            ignore
                Accept invalid residues without changes.

            fail
                Raise an exception on the first invalid residue.

            remove
                Discard any sequence that contains invalid residues.

            convert-all
                Apply the standard conversion table
                ``B->N, U->C, X->G, Z->Q, '*'->'', '-'->''``
                and raise an exception if any residues remain invalid
                afterwards.

            convert-res
                Same as ``convert-all`` but keeps the alignment gap
                character ``'-'`` untouched.

            convert-all-ignore
                Same as ``convert-all`` but silently keeps any residues
                that remain invalid after conversion.

            convert-res-ignore
                Same as ``convert-res`` but silently keeps any residues
                that remain invalid after conversion.

            convert-all-remove
                Same as ``convert-all`` but removes any sequence that
                still contains invalid residues after conversion.

            convert-res-remove
                Same as ``convert-res`` but removes any sequence that
                still contains invalid residues after conversion.

    --number-lines <int>                     (default: 60)
        Number of residues per line in the output FASTA file. Must be
        at least 5.

    --shortest-seq <int>                     (default: none)
        Minimum sequence length to include. Sequences shorter than or
        equal to this length are discarded.

    --longest-seq <int>                      (default: none)
        Maximum sequence length to include. Sequences longer than or
        equal to this length are discarded. If both ``--longest-seq``
        and ``--shortest-seq`` are given, ``--longest-seq`` must be
        larger.

    --random-subsample <int>                 (default: none)
        Randomly sub-sample this many sequences from the final set.
        Useful for generating small test FASTA files from large inputs.
        If the input contains fewer sequences than requested, all
        sequences are returned.

    --print-statistics
        Print length statistics (count, 25th / 50th / 75th percentile,
        longest, shortest) for the final set of sequences.

    --no-outputfile
        Do not write an output FASTA file. Useful together with
        ``--print-statistics`` for pure summary runs.

    --remove-comma-from-header
        Replace ``,`` with ``;`` in every FASTA header on read. Useful
        when downstream tools parse FASTA headers as CSV fields.

    --silent
        Suppress all ``[INFO]`` output to stdout.

    --version
        Print the installed **protfasta** version and exit.


Examples
.........

Clean up a FASTA file by removing duplicate records and converting
non-standard residues, writing the result to ``clean.fasta``::

    pfasta --duplicate-record remove --invalid-sequence convert-all \
           -o clean.fasta input.fasta

Filter sequences between 50 and 500 residues and randomly keep 1000
of them::

    pfasta --shortest-seq 50 --longest-seq 500 \
           --random-subsample 1000 \
           -o subset.fasta input.fasta

Just print length statistics, without writing a file::

    pfasta --print-statistics --no-outputfile input.fasta


.. toctree::
   :maxdepth: 2
   :caption: Contents:
