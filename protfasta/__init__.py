"""
protfasta - A simple but robuts FASTA parser explicitly for protein sequences.


.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other. 

"""


# Add imports here
from .protfasta import read_fasta
from ._configs import STANDARD_AAS, STANDARD_CONVERSION
from . import io as _io
from .io import  write_fasta_file 
import os

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


_ROOT = os.path.abspath(os.path.dirname(__file__))
def _get_data(path):    
    """
    This function is used for getting the absolute path for loading
    test data. Pay no attention to the code behind the curtain.
    
    """
    return os.path.join(_ROOT, 'data', path)
