# Package data

`test_data/` holds the small FASTA files the test suite reads: a clean UniProt-style set, files with duplicate records and duplicate sequences, files with convertible and unconvertible non-standard residues, and a few short alignments. They are installed with the package (see `graft protfasta` in `MANIFEST.in`) so that the suite can be run against an installed copy with `pytest --pyargs protfasta.tests`, which is what tox and the CI matrix do. Tests locate the directory through `protfasta._get_data('test_data')` rather than a path relative to the source tree.

Anything the tests write goes to pytest-managed temporary directories, so nothing under here is modified by a test run.
