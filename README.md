protfasta
==============================
[//]: # (Badges)
[![PyPI version](https://img.shields.io/pypi/v/protfasta.svg)](https://pypi.org/project/protfasta/)
[![Python versions](https://img.shields.io/pypi/pyversions/protfasta.svg)](https://pypi.org/project/protfasta/)
[![Documentation Status](https://readthedocs.org/projects/protfasta/badge/?version=latest)](https://protfasta.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/github/license/holehouse-lab/protfasta.svg)](https://github.com/holehouse-lab/protfasta/blob/master/LICENSE)
[![Last commit](https://img.shields.io/github/last-commit/holehouse-lab/protfasta.svg)](https://github.com/holehouse-lab/protfasta/commits/master)
[![Open issues](https://img.shields.io/github/issues/holehouse-lab/protfasta.svg)](https://github.com/holehouse-lab/protfasta/issues)
[![Stars](https://img.shields.io/github/stars/holehouse-lab/protfasta.svg?style=flat)](https://github.com/holehouse-lab/protfasta/stargazers)
[![Downloads](https://img.shields.io/pypi/dm/protfasta.svg)](https://pypi.org/project/protfasta/)



## Release 0.1.26 (October 2026)

## Overview
protfasta - a robust parser for protein-based FASTA files.

## Documentation

For all documentation see [https://protfasta.readthedocs.io/en/latest/](https://protfasta.readthedocs.io/en/latest/).

For code see [https://github.com/holehouse-lab/protfasta](https://github.com/holehouse-lab/protfasta).

If you want to change the package rather than use it, [ARCHITECTURE.md](ARCHITECTURE.md) walks through how a `read_fasta` and a `write_fasta` call move through the modules.

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

For files that are too large to fit in memory, `read_fasta_stream()` applies the same sanitization but yields one record at a time:

	import protfasta
	
	for header, sequence in protfasta.read_fasta_stream('huge.fasta'):
	    ...

## Running tests

To run tests on your environment, clone the source code and from protfasta root 

```bash
cd protfasta/tests
pytest --verbose
```

To run the test suite over all supported Python versions we use [tox](https://tox.wiki/en/4.58.0/);  from the protfasta source root directory simply run

```bash
tox	
```

And the tests should run across Python envs 3.9 to 3.15 (experimental).



## Errors and help

For bug reports or errors please raise an issue on this github repository (see the [Issues](https://github.com/holehouse-lab/protfasta/issues) tab at the top).

## Changelog

The full release history lives in [CHANGELOG.md](CHANGELOG.md).

## Copyright

Copyright (c) 2020-2026, Alex Holehouse  - [Holehouse lab](http://holehouse.wustl.edu/). `protfasta` is released under the MIT license. The codebase is well structured and relatively simple, lending it to feature addition. We welcome pull-requests assuming contributed code maintains an appropriate level of clarity and robustness. 


#### Acknowledgements

Many of the software-engineering tools and approaches used in the development of `protfasta` come from resources developed by the [Molecular Sciences Software Institute](https://molssi.org/).
