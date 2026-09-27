Installation
===============

**protfasta** requires Python 3.9 or higher and has no dependencies beyond
the Python standard library. It has been tested on Linux and macOS. It should also work on Windows but we haven't tested it there yet.

**protfasta** can be downloaded and installed directly from PyPI using **pip**:

.. code-block:: none

   pip install protfasta


If this has worked, the **pfasta** tool should be available from the command-line

.. code-block:: none

   pfasta --help


And you're done. This also means you can now ``import`` and use **protfasta** in your Python workflow.


Installing from source
.......................

To install the current development version instead, clone the repository and install it with **pip**:

.. code-block:: none

   git clone https://github.com/holehouse-lab/protfasta.git
   cd protfasta
   pip install .


Running the tests
..................

The test suite is installed along with the package, so once **protfasta** is installed (from PyPI or from source) it can be run from anywhere:

.. code-block:: none

   pip install pytest
   pytest --pyargs protfasta.tests

From a clone of the repository, running ``tox`` in the repository root runs the same suite across every supported Python version.
