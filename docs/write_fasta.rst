write_fasta
=================

``write_fasta`` writes a set of sequences to a FASTA file. It accepts
sequence data as either:

    *  a dictionary ``{header: sequence, ...}``, or
    *  a list of ``[header, sequence]`` pairs.

It is also possible to have :func:`protfasta.read_fasta` write its
sanitized result directly to disk via the ``output_filename`` keyword,
which simply calls ``write_fasta`` internally.


Keyword arguments
...................

    *  ``filename`` - destination path, either a string or a
       :class:`pathlib.Path`. Conventionally ends with ``.fasta`` or
       ``.fa`` but this is not enforced.

    *  ``linelength`` (default ``60``) - maximum residues per line.
       Values from ``1`` to ``4`` are clamped to ``5``. Set to ``0`` (or
       any negative value), ``None`` or ``False`` to write each sequence
       on a single line. A numerical string such as ``'60'`` is accepted
       and cast.

    *  ``append_to_fasta`` (default ``False``) - when ``True``, new
       entries are appended to an existing file rather than
       overwriting it. Must be a boolean.


Error handling
...............

``write_fasta`` raises a ``ProtfastaException`` rather than letting a
lower-level error escape, so a single ``except`` clause is enough:

    *  ``fasta_data`` is neither a dictionary nor a list.
    *  An element of a list is not a two-item ``[header, sequence]``
       pair (a two-item tuple is accepted as well).
    *  A header or sequence is not a string.
    *  A header or sequence contains a line break, which would be read
       back as a record boundary and silently corrupt the file.
    *  A sequence is empty - empty records are never silently written.
    *  ``linelength`` cannot be interpreted as an integer. Note that a
       numerical string (``'60'``) *is* accepted and cast.
    *  ``append_to_fasta`` is not a boolean.
    *  ``filename`` cannot be opened for writing (for example, its
       directory does not exist).

All of the data checks run **before** the file is opened, so a bad entry
never leaves a truncated file behind and never partially appends to an
existing one.


Encoding
.........

Output is written as UTF-8. Headers that came from a non-UTF-8 file via
:func:`protfasta.read_fasta` (see its Encoding section) are written back
out with their original bytes, so a read/write round trip is lossless.


Performance notes
..................

``write_fasta`` assembles each record into a single string and emits it
with one write call into a 1 MiB buffer, which makes it suitable for
very large outputs (tens of millions of sequences and beyond).


For usage examples see the :doc:`examples` page.


Documentation
...............


.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: protfasta
   :noindex:

.. autofunction:: write_fasta
