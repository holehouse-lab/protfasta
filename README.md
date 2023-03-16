protfasta
==============================
[//]: # (Badges)
[![Build Status](https://travis-ci.org/holehouse-lab/protfasta.svg?branch=master)](https://travis-ci.org/holehouse-lab/protfasta)

## Release 0.1.12 (March 2023)

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
* **0.1.12** (March 2023) - integrated in check_header_parser flag via pull request from the amazing [Friedlab](https://friedlab.com/) !
* Added in `append_to_fasta` flag so you can append to an existing FASTA file (thanks Ryan!)

* **0.1.11** (Sept 17th 2022) - re-wrote code for checking duplicate sequence to make it O(1) instead of O(n) for number of sequences (:-/) and added convert-remove option for invalid_sequences

* **0.1.9** (Sept 12th 2021) - added in robustness for whitespace in sequence files, which, bizarrely, was not present (i.e. added as an invalid residue type but can now be converted).


## Copyright

Copyright (c) 2020-2021, Alex Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/). `protfasta` is released under the MIT license. The codebase is well structured and relatively simple, lending it to feature addition. We welcome pull-requests assuming contributed code maintains an appropriate level of clarity and robustness. 


#### Acknowledgements
 
Many of the software-engineering tools and approaches used in the development of `protfasta` come from resources developed by the [Molecular Sciences Software Institute](https://molssi.org/).
