"""
Custom exception classes for protfasta.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""


class ProtfastaException(Exception):
    """Exception raised for protfasta-specific errors.

    Inherits from :class:`Exception` and is used throughout the package
    to signal invalid input, duplicate records, unrecognised keyword
    arguments, and other domain-specific error conditions.
    """
