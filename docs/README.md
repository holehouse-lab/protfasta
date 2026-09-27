# Compiling protfasta's Documentation

The docs for this project are built with [Sphinx](https://www.sphinx-doc.org/) using the Read the Docs theme, and are published at [protfasta.readthedocs.io](https://protfasta.readthedocs.io/). To build them locally, install the documentation requirements (which pull in Sphinx itself) from this directory:

```bash
pip install -r requirements.txt
```

Then either use the `Makefile`:

```bash
make html
```

or call Sphinx directly, which is what the Read the Docs build does:

```bash
python -m sphinx -b html . _build/html
```

Adding `-W` turns warnings into errors, which is a useful check before a release. The compiled pages end up in `_build/html`; open `_build/html/index.html` in a browser to view them. The API sections are generated from the docstrings in the `protfasta` package, and `conf.py` puts the repository root on the path, so the build documents the checked-out code rather than any installed copy.
