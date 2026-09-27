#!/usr/bin/env python
"""Generate a synthetic protein FASTA file with a given number of records.

The point of this script is to produce large, realistic-looking input for the
protfasta benchmarks quickly and reproducibly. Sequence lengths are drawn
from a log-normal distribution (median ~350 residues, clipped to 30-5000),
which is a reasonable stand-in for a proteome; residues are drawn uniformly
from the 20 standard amino acids and headers mimic the UniProt format.
Sequences are wrapped at 60 residues per line, as UniProt does.

Generation is vectorised with numpy so that a 10-million-record file (~5 GB)
takes on the order of a minute rather than the better part of an hour.

Usage
-----
    python generate_fasta.py OUTPUT.fasta N_RECORDS [--seed 0] [--linelength 60]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

AAS = b'ACDEFGHIKLMNPQRSTVWY'

# log-normal length parameters (in residues)
_MEDIAN_LENGTH = 350
_SIGMA = 0.6
_MIN_LENGTH = 30
_MAX_LENGTH = 5000


def generate_fasta(
    filename: str,
    n_records: int,
    seed: int = 0,
    linelength: int = 60,
    chunk: int = 200_000,
) -> int:
    """Write a synthetic FASTA file with *n_records* records.

    Parameters
    ----------
    filename : str
        Output path (overwritten if it exists).

    n_records : int
        Number of records to write.

    seed : int, optional
        Seed for the random generator; the same seed always produces the
        same file. Default ``0``.

    linelength : int, optional
        Residues per sequence line. ``0`` writes each sequence on a single
        line. Default ``60``.

    chunk : int, optional
        Records generated per numpy batch. Only affects peak memory of the
        generator, not the output. Default ``200000``.

    Returns
    -------
    int
        Total number of residues written.
    """
    rng = np.random.default_rng(seed)
    table = np.frombuffer(AAS, dtype=np.uint8)

    total_residues = 0
    written = 0
    with open(filename, 'wb', buffering=8 << 20) as fh:
        while written < n_records:
            m = min(chunk, n_records - written)

            lengths = rng.lognormal(mean=np.log(_MEDIAN_LENGTH), sigma=_SIGMA, size=m)
            lengths = np.clip(lengths.astype(np.int64), _MIN_LENGTH, _MAX_LENGTH)
            n_res = int(lengths.sum())

            # one contiguous buffer of residues for the whole batch. The
            # indices only need to span 0-19, so drawing them as uint8 rather
            # than the default int64 cuts the generator's peak memory by
            # ~600 MB per batch of 200,000 records.
            residues = np.take(table, rng.integers(0, len(AAS), size=n_res, dtype=np.uint8))
            buf = residues.tobytes()
            offsets = np.concatenate(([0], np.cumsum(lengths)))

            parts: list[bytes] = []
            for i in range(m):
                idx = written + i
                seq = buf[offsets[i]:offsets[i + 1]]
                parts.append(b'>sp|P%08d|PROT%d_HUMAN Synthetic protein %d OS=Homo sapiens OX=9606 GN=G%d PE=1 SV=1\n'
                             % (idx, idx, idx, idx))
                if linelength:
                    parts.append(b'\n'.join([seq[j:j + linelength] for j in range(0, len(seq), linelength)]))
                else:
                    parts.append(seq)
                parts.append(b'\n')

            fh.write(b''.join(parts))
            written += m
            total_residues += n_res

    return total_residues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output', help='output FASTA filename')
    parser.add_argument('n_records', type=int, help='number of records to generate')
    parser.add_argument('--seed', type=int, default=0, help='random seed (default 0)')
    parser.add_argument('--linelength', type=int, default=60, help='residues per line, 0 for single-line sequences (default 60)')
    args = parser.parse_args()

    t0 = time.perf_counter()
    n_res = generate_fasta(args.output, args.n_records, seed=args.seed, linelength=args.linelength)
    dt = time.perf_counter() - t0
    size_mb = os.path.getsize(args.output) / 1e6
    print('Wrote %s: %d records, %d residues, %.1f MB in %.1f s' % (args.output, args.n_records, n_res, size_mb, dt),
          file=sys.stderr)


if __name__ == '__main__':
    main()
