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

from typing import Mapping, Union

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

# Internal frozensets used for O(1) membership tests (public list versions
# are kept above for backwards compatibility).
_STANDARD_AAS_SET: frozenset[str] = frozenset(STANDARD_AAS)
_STANDARD_AAS_WITH_GAP_SET: frozenset[str] = frozenset(STANDARD_AAS_WITH_GAP)

# Pre-built translation tables for fast single-pass sequence cleaning via
# ``str.translate``.  Built once at import time to avoid per-call rebuild.
# (Keyed by code point, which is exactly what str.maketrans would produce
# for single-character keys.)
_TRANSLATE_STANDARD: Mapping[int, Union[str, int, None]] = {ord(k): v for k, v in STANDARD_CONVERSION.items()}
_TRANSLATE_WITH_GAP: Mapping[int, Union[str, int, None]] = {ord(k): v for k, v in STANDARD_CONVERSION_WITH_GAP.items()}

# Tables used to *validate* a sequence in a single C-level pass. Deleting
# every valid residue from a sequence leaves behind exactly the invalid
# characters (in order), so an empty result means the sequence is clean.
#
# The ``bytes`` variants are used for the (overwhelmingly common) pure-ASCII
# case, where ``bytes.translate(None, valid)`` is a tight table-driven loop
# that runs several times faster than the ``str`` equivalent; the ``str``
# tables handle anything containing non-ASCII characters.
_VALID_BYTES_STANDARD: bytes = ''.join(STANDARD_AAS).encode('ascii')
_VALID_BYTES_WITH_GAP: bytes = ''.join(STANDARD_AAS_WITH_GAP).encode('ascii')
_DELETE_VALID_STANDARD = str.maketrans('', '', ''.join(STANDARD_AAS))
_DELETE_VALID_WITH_GAP = str.maketrans('', '', ''.join(STANDARD_AAS_WITH_GAP))
