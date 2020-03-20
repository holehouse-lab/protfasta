"""
protfasta
A hyper robuts FASTA parser explicitly for protein sequences
"""

# Add imports here
from .protfasta import read_fasta
import os

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)
