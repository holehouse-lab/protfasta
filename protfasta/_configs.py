"""
Internal configuration constants for protfasta.

Defines the standard amino-acid alphabet, the extended alphabet that
includes gap characters, and the default character-conversion tables
used when cleaning non-standard residues from protein sequences.

.............................................................................
protfasta was developed by the Holehouse lab
     Original release March 2020

Question/comments/concerns? Raise an issue on github:
https://github.com/holehouse-lab/protfasta

Licensed under the MIT license.

Be kind to each other.

"""

from __future__ import annotations

STANDARD_CONVERSION: dict[str, str] = {
    'B': 'N',
    'U': 'C',
    'X': 'G',
    'Z': 'Q',
    '*': '',
    '-': '',
    ' ': '',
}
"""Default mapping of non-standard characters to standard replacements."""

STANDARD_CONVERSION_WITH_GAP: dict[str, str] = {
    'B': 'N',
    'U': 'C',
    'X': 'G',
    'Z': 'Q',
    ' ': '',
    '*': '',
}
"""Like :data:`STANDARD_CONVERSION` but preserves dashes as gap characters."""

STANDARD_AAS: list[str] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
]
"""The 20 standard amino-acid one-letter codes."""

STANDARD_AAS_WITH_GAP: list[str] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-',
]
"""The 20 standard amino acids plus the dash gap character."""
