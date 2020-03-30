.. protfasta documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

protfasta - the simple FASTA parser for proteins
=========================================================

**protfasta** is a simple, robust parser for working with FASTA files. It is pure python and has no external package dependencies other than Python language modules.

It contains two distinct components:

    1. A Python API for reading and writing FASTA files, which includes a collection of santization functions. This makes it easy to write code that reads/writes FASTA files.

    2. A command line tool (``pfasta``) that allows manipulation of FASTA files directly from the command line.

This documentation provides an overeview of both components.

Why did you build protfasta? Don't you have better things to do with your time?
..............
This is a reasonable question...

Working with protein-based FASTA files is at the heart of a lot of what `the Holehouse lab <http://holehouse.wustl.edu/>`_ does. We had previously used a few different existing parsers but had found limitations with respect to certain features. Part of this came from the fact that many FASTA parsers can work with nucleotide or protein data. Given our bread and butter is protein sequences, we decided to build a parser explicitly for working with proteins. We also wanted the ability to deal with FASTA files with duplicate entries. Not necessarily because this is 'good', but for processing reasons being able to deal with this in-code is easier than having to sanitize ahead of time.

We built ``pfasta`` as a compact tool for working with FASTA files at the command-line level. In particular, the ability to filter FASTA files by sequence length, correct/remove sequence with invalid amino acids, and do various other things lends ``pfasta`` as a useful first tool in our informatics pipelines.


Will protfasta work with nucleotide-based FASTA files?
.......................
In principle yes, but none of our testing suites are set up to rigerously explore this. However, there's no reason it shouldn't, although it may be less efficient that some other tools such as the excellent `pyfaidx <https://pypi.org/project/pyfaidx/>`_.


Bugs and help
..............

If you find any bugs or have feature requests please raise an issue on our `Github page <https://github.com/holehouse-lab/protfasta/>`_. **protfasta** uses a continous integration suite for the main package, and **pfasta** has a set of local tests that are run upon updates.


How to cite **protfasta**
.........................
For now please just cite the `Github repository <https://github.com/holehouse-lab/protfasta/>`_ (including the date accessed). We are planning on putting out a short biorxiv paper (not submitting to a journal) for a DOI-ed reference, at which point this documentation will be updated accordingly.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   pfasta
   examples
   read_fasta
   write_fasta





