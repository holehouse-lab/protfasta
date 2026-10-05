pfasta
=================

**pfasta** is a command-line tool for working with FASTA files to filter and
sanitize them based on various criteria. It is installed automatically with
**protfasta** and can be invoked as ``pfasta`` from the command line (or as
``python -m protfasta.scripts.pfasta``).

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

Under the hood **pfasta** reads the file with :func:`protfasta.read_fasta`
and writes it with :func:`protfasta.write_fasta`, so everything those pages
say about parsing, sanitization and the output format applies here too.

Usage
.........

.. code-block:: none

    pfasta <flags> filename.fasta


Command-line options
.....................

.. code-block:: none

    filename
        Positional argument: path to the input FASTA file.

    -h, --help
        Print a summary of the options and exit.

    -o <output filename>                    (default: output.fasta)
        Output FASTA file. It is created (or overwritten) in the current
        directory if no path is given.

    --non-unique-header
        If set, multiple FASTA records are allowed to share the same header.
        By default a duplicate header is a fatal error.

    --duplicate-record {ignore,fail,remove} (default: fail)
        How to deal with duplicate records (same header AND same sequence):
            fail   - report the first duplicate record and exit
            ignore - keep all duplicate records
            remove - keep only the first occurrence
        A duplicate record can only exist if headers are allowed to repeat,
        so 'ignore' requires --non-unique-header, and 'fail' and 'remove'
        only have anything to act on when it is set.

    --duplicate-sequence {ignore,fail,remove} (default: ignore)
        How to deal with duplicate sequences (same sequence, any header):
            fail   - report the first duplicate sequence and exit
            ignore - keep all duplicate sequences
            remove - keep only the first occurrence of each sequence

    --empty-sequence {fail,remove}            (default: fail)
        How to deal with a header that has no sequence after it:
            fail   - report the header and exit
            remove - drop the record
        ('ignore' is not offered, because a record with no sequence
        cannot be written to the output FASTA file.)

    --invalid-sequence <mode>                (default: fail)
        How to deal with non-standard amino-acid characters. Available
        modes:

            ignore
                Accept invalid residues without changes.

            fail
                Report the first invalid residue and exit.

            remove
                Discard any sequence that contains invalid residues.

            convert-all
                Apply the standard conversion table
                B->N, U->C, X->G, Z->Q, '*'->'', ' '->'', '-'->''
                and exit with an error if any residues remain invalid
                afterwards.

            convert-res
                Same as convert-all but the file is read as an
                alignment: the gap character '-' is preserved and
                treated as a valid residue rather than converted.

            convert-all-ignore
                Same as convert-all but silently keeps any residues
                that remain invalid after conversion.

            convert-res-ignore
                Same as convert-res but silently keeps any residues
                that remain invalid after conversion.

            convert-all-remove
                Same as convert-all but removes any sequence that
                still contains invalid residues after conversion.

            convert-res-remove
                Same as convert-res but removes any sequence that
                still contains invalid residues after conversion.

    --number-lines <int>                     (default: 60)
        Number of residues per line in the output FASTA file. Must be
        at least 5.

    --shortest-seq <int>                     (default: none)
        Minimum sequence length to include. Sequences shorter than or
        equal to this length are discarded.

    --longest-seq <int>                      (default: none)
        Maximum sequence length to include. Sequences longer than or
        equal to this length are discarded. If both --longest-seq
        and --shortest-seq are given, --longest-seq must be larger.

    --random-subsample <int>                 (default: none)
        Randomly sub-sample this many sequences from the final set
        (after the length filters). Useful for generating small test
        FASTA files from large inputs. The sub-sample is written in
        random order rather than file order; if the input contains
        fewer sequences than requested, all of them are written, again
        in random order.

    --print-statistics
        Print length statistics (count, 25th / 50th / 75th percentile,
        longest, shortest) for the final set of sequences. Nothing is
        printed if --silent is also given.

    --no-outputfile
        Do not write an output FASTA file. Useful together with
        --print-statistics for pure summary runs.

    --remove-comma-from-header
        Replace ',' with ';' in every FASTA header on read. Useful
        when downstream tools parse FASTA headers as CSV fields.

    --silent
        Print nothing to stdout except a fatal error: no banner, no
        [INFO] progress lines and no --print-statistics output.

    --version
        Print the installed protfasta version and exit.


Examples
.........

Clean up a FASTA file by removing duplicate records and converting
non-standard residues, writing the result to ``clean.fasta``. Duplicate
records share a header, so ``--non-unique-header`` is needed for them to
reach the removal step at all (without it the first repeated header is a
fatal error)::

    pfasta --non-unique-header --duplicate-record remove \
           --invalid-sequence convert-all \
           -o clean.fasta input.fasta

Keep sequences that are longer than 50 and shorter than 500 residues and
randomly keep 1000 of them::

    pfasta --shortest-seq 50 --longest-seq 500 \
           --random-subsample 1000 \
           -o subset.fasta input.fasta

Just print length statistics, without writing a file::

    pfasta --print-statistics --no-outputfile input.fasta


Exit status
............

If the input cannot be read or sanitized under the requested options (a
duplicate header, a duplicate record under ``--duplicate-record fail``, an
invalid residue under ``--invalid-sequence fail``, an output path that
cannot be created, and so on), or an option has an invalid value (an
unknown ``--invalid-sequence`` mode, a non-numeric length), **pfasta**
prints a single ``[FATAL ERROR]`` line describing the problem and exits with
status ``1``. No output file is written in that case.

An unrecognised flag, or a missing input filename, is reported by the
argument parser on stderr with a short usage summary, and **pfasta** exits
with status ``2``.

If the length filters or ``--random-subsample`` leave no sequences at
all, **pfasta** reports ``0 sequences remain after filtering`` and exits
with status ``0``, again without writing an output file.
