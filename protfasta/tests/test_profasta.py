"""
Unit and regression test for the protfasta package.
"""

# Import package, test suite, and other packages as needed
import protfasta
import pytest
import sys

def test_protfasta_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "protfasta" in sys.modules
