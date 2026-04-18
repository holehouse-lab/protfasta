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

    *  ``filename`` - destination path. Conventionally ends with
       ``.fasta`` or ``.fa`` but this is not enforced.

    *  ``linelength`` (default ``60``) - maximum residues per line.
       Values below ``5`` are clamped to ``5``. Set to ``0``, ``None``
       or ``False`` to write each sequence on a single line.

    *  ``append_to_fasta`` (default ``False``) - when ``True``, new
       entries are appended to an existing file rather than
       overwriting it.


Performance notes
..................

``write_fasta`` writes output in line-length chunks with a 1 MiB write
buffer, which makes it suitable for very large outputs (tens of
millions of sequences and beyond). An empty sequence raises a
``ProtfastaException`` rather than being silently written out.


For usage examples see the :doc:`examples` page.


Documentation
...............


.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: protfasta
   :noindex:

.. autofunction:: write_fasta
