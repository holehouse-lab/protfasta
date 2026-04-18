read_fasta
=================

``read_fasta`` is the primary entry point to **protfasta**. It reads a
FASTA file, optionally sanitizes its contents, and returns the records
either as a dictionary (``header -> sequence``) or as a list of
``[header, sequence]`` pairs.

At its simplest::

    import protfasta
    sequences = protfasta.read_fasta('proteins.fasta')

Many optional keyword arguments customize the behaviour - including
duplicate handling, invalid-residue correction, alignment-gap support,
custom header parsing, and automatic writing of the sanitized output.


What ``read_fasta`` can do
...........................

    *  Ignore, remove, convert, or fail on sequences containing
       non-standard amino-acid characters
       (``B``, ``U``, ``X``, ``Z``, ``*``, ``-``, ``' '``, ...).
    *  Ignore, remove, or fail on duplicate FASTA records
       (same header **and** same sequence).
    *  Ignore, remove, or fail on duplicate sequences (same sequence,
       different headers).
    *  Preserve alignment gap characters (``-``) when
       ``alignment=True``.
    *  Apply a caller-supplied ``header_parser`` function to every
       raw header (useful for extracting accession IDs).
    *  Optionally write the sanitized result to a new FASTA file via
       ``output_filename``.
    *  Override the built-in invalid-character conversion table via
       a custom ``correction_dictionary``.


Processing pipeline
.....................

Sanitization happens in a fixed order:

    1. File is streamed from disk, headers are parsed with
       ``header_parser`` (if provided), and header uniqueness is
       checked (when ``expect_unique_header=True``).
    2. Duplicate **records** are processed according to
       ``duplicate_record_action``.
    3. Duplicate **sequences** are processed according to
       ``duplicate_sequence_action``.
    4. **Invalid residues** are processed according to
       ``invalid_sequence_action``.
    5. The sanitized set is optionally written to
       ``output_filename``.
    6. The result is returned as a ``dict`` (default) or a ``list``
       of ``[header, sequence]`` pairs (when ``return_list=True``).

Incompatible option combinations (for example,
``expect_unique_header=True`` together with
``duplicate_record_action='ignore'``) are caught before the file is
read.


Default conversion table
..........................

When ``invalid_sequence_action`` includes conversion and no custom
``correction_dictionary`` is supplied, these replacements are applied:

    *  ``B`` -> ``N``
    *  ``U`` -> ``C``
    *  ``X`` -> ``G``
    *  ``Z`` -> ``Q``
    *  ``*`` -> ``''`` (removed)
    *  ``-`` -> ``''`` (removed; preserved if ``alignment=True``)
    *  ``' '`` -> ``''`` (whitespace removed)


Large files
............

For files that do not fit comfortably in memory, consider using the
streaming parser :func:`protfasta.iter_fasta` instead. ``read_fasta``
itself streams the file from disk (so it will not load the entire file
as a single string), but it still builds an in-memory data structure
of all records; ``iter_fasta`` avoids that.


For usage examples see the :doc:`examples` page. Full API
documentation is shown below.


Documentation
...............

.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. automodule:: protfasta

.. autofunction:: read_fasta
