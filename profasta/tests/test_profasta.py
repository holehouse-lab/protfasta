"""
Unit and regression test for the profasta package.
"""

# Import package, test suite, and other packages as needed
import profasta
import pytest
import sys

def test_profasta_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "profasta" in sys.modules
