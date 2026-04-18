protfasta
==============================
[//]: # (Badges)
[![Build Status](https://travis-ci.org/holehouse-lab/protfasta.svg?branch=master)](https://travis-ci.org/holehouse-lab/protfasta)



## Release 0.1.15 (October 2024)

## Overview
protfasta - a robust parser for protein-based FASTA files.

## Documentation

For all documentation see [https://protfasta.readthedocs.io/en/latest/](https://protfasta.readthedocs.io/en/latest/).

For code see [https://github.com/holehouse-lab/protfasta](https://github.com/holehouse-lab/protfasta).

## Installation

`protfasta` has been tested on Linux and macOS. It should also work on Windows but we haven't tested it there yet. 

`protfasta` can be downloaded and installed directly from PyPI using **pip**:

    pip install protfasta

If this has worked, the `pfasta` command-line tool should be available from the command-line

    pfasta --help

And you're done. This also means you can now ``import`` and use **protfasta** in your Python workflow. 

## Simple example

	import protfasta
	
	# sequences is now a dictionary where keys are FASTA headers and values are sequences.
	sequences = protfasta.read_fasta('inputfile.fasta')


## Errors and help
For bug reports or errors please raise an issue on this github repository (see the [Issues](https://github.com/holehouse-lab/protfasta/issues) tab at the top).

## Changelog
* **0.1.16** (April 2026) - Performance overhaul for large FASTA files (hundreds of millions of sequences).
	* `read_fasta` now streams the input file instead of reading it entirely into memory with `readlines()` — peak memory is now O(single record) rather than O(file size).
	* New `protfasta.iter_fasta(filename, header_parser=None)` generator for memory-bounded streaming access to `(header, sequence)` pairs from files that don't fit in RAM.
	* Core parser rewritten to use list-of-parts + `''.join()` instead of quadratic string concatenation, and to skip header-uniqueness tracking entirely when `expect_unique_header=False`.
	* `convert_to_valid` / invalid-residue handling now uses a pre-built `str.translate` table (single C-level pass) instead of a chained `str.replace` loop — typically 5–20× faster.
	* `check_sequence_is_valid` now uses frozenset membership instead of list-based `in` scans.
	* Duplicate-detection utilities (`fail_on_duplicates`, `remove_duplicates`, `fail_on_duplicate_sequences`, `remove_duplicate_sequences`) now store 16-byte blake2b digests instead of full sequences in their lookup structures, dramatically reducing peak memory for files with long sequences.
	* `write_fasta` replaced its per-residue `fh.write()` loop with chunked slice writes and opens the output with a 1 MiB buffer — roughly two orders of magnitude faster on large files.
	* All existing behavior and the full test suite (239 tests) are preserved.

* **0.1.14**  and **0.1.15** (October 2024) - Re-wrote build chain and versioning to use `pyproject.toml` and [versioningit](https://pypi.org/project/versioningit/). protfasta should now support Python beyond 3.12. About bloody time. 
	* Added `--version` flag to pfasta
	* Messed around a bit with tags to ensure we had a tagged version compatible with them. 

* **0.1.13** (January 2023) - Added upper limit of Python 3.11 to accomodate clash between versioneer and Python 3.12. Ultimately we'll move to versioningit for release versioning (as we have done internally) but need to make sure we have a robust protocol for this switch and then do this for ALL tools....

* **0.1.12** (March 2023) - integrated in check_header_parser flag via pull request from the amazing [Friedlab](https://friedlab.com/) !
* Added in `append_to_fasta` flag so you can append to an existing FASTA file (thanks Ryan!)

* **0.1.11** (Sept 17th 2022) - re-wrote code for checking duplicate sequence to make it O(1) instead of O(n) for number of sequences (:-/) and added convert-remove option for invalid_sequences

* **0.1.9** (Sept 12th 2021) - added in robustness for whitespace in sequence files, which, bizarrely, was not present (i.e. added as an invalid residue type but can now be converted).


## Copyright

Copyright (c) 2020-2026, Alex Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/). `protfasta` is released under the MIT license. The codebase is well structured and relatively simple, lending it to feature addition. We welcome pull-requests assuming contributed code maintains an appropriate level of clarity and robustness. 


#### Acknowledgements

Many of the software-engineering tools and approaches used in the development of `protfasta` come from resources developed by the [Molecular Sciences Software Institute](https://molssi.org/).
