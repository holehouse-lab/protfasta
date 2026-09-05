#!/usr/bin/env python
"""Benchmark protfasta's parsing speed and memory footprint against file size.

For each requested number of records a synthetic FASTA file is generated
(once, and cached in ``--data-dir``) with :mod:`generate_fasta`, and then a
set of read modes is timed:

* ``read_fasta``           - ``protfasta.read_fasta`` with default options
                             (unique headers enforced, invalid residues fail).
* ``read_fasta-nocheck``   - ``read_fasta`` with every check disabled, i.e.
                             raw parse cost plus building the result.
* ``stream``               - ``protfasta.read_fasta_stream`` with default
                             options, fully consumed. This is the flat-memory
                             configuration.
* ``stream-unique``        - ``read_fasta_stream(expect_unique_header=True)``,
                             which keeps a set of headers (O(records) memory).
* ``stream-dedup``         - ``read_fasta_stream(duplicate_sequence_action=
                             'remove')``, which keeps a set of 16-byte
                             sequence digests (O(records) memory).

Every (size, mode) measurement runs in a **fresh subprocess** so that peak
memory is not contaminated by earlier runs. Before the modes for a given
size are timed the input file is read once, sequentially, to bring it into
the page cache; without this the first mode measured after generating (or
after a reboot) pays for disk I/O that the others do not. Memory is reported as the
subprocess's peak resident set size (``ru_maxrss``) minus the peak RSS of a
subprocess that only imports protfasta, so the figure is the memory the
parse itself needed. For ``read_fasta`` that includes the returned data
structure - the whole point of the streaming modes is that it does not.

Timing covers only the read call, not interpreter start-up or imports.
With ``--repeats N`` the fastest time and largest memory over N runs are
kept.

Usage
-----
    python benchmark_protfasta.py [--sizes 10000 100000 1000000 10000000]
                                  [--modes read_fasta stream ...]
                                  [--repeats 1] [--data-dir DIR] [--results-dir DIR]

Results are printed as a Markdown table and also written as JSON and
Markdown to ``--results-dir``.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import resource
import subprocess
import sys
import time
from typing import Any, Callable, Optional, cast

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)

# Make sure the worker subprocesses import the checkout, not an installed copy.
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, HERE)

MODES: list[str] = ['read_fasta', 'read_fasta-nocheck', 'stream', 'stream-unique', 'stream-dedup']
DEFAULT_SIZES: list[int] = [10_000, 100_000, 1_000_000, 10_000_000]


# ---------------------------------------------------------------------------
# worker side: runs inside a fresh subprocess
# ---------------------------------------------------------------------------
def _peak_rss_bytes() -> int:
    """Peak resident set size of this process, in bytes, on any platform."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes, Linux reports kilobytes
    return rss if sys.platform == 'darwin' else rss * 1024


def _run_mode(mode: str, filename: str) -> dict[str, Any]:
    """Execute one read mode on *filename* and return timing/memory data."""
    import protfasta

    n_records = 0
    n_residues = 0

    # For the read_fasta modes the clock stops as soon as the call returns;
    # the record and residue counts are taken afterwards so that the extra
    # pass over the result is not charged to the parse. The streaming modes
    # count inside the consuming loop, which is the work being measured.
    t0 = time.perf_counter()
    if mode == 'import-only':
        pass

    elif mode == 'read_fasta':
        # (return_list is False, so the result is a dict - cast for the type checker)
        data = cast(dict, protfasta.read_fasta(filename))
        seconds = time.perf_counter() - t0
        n_records = len(data)
        n_residues = sum(map(len, data.values()))

    elif mode == 'read_fasta-nocheck':
        data = cast(dict, protfasta.read_fasta(filename,
                                               expect_unique_header=False,
                                               duplicate_record_action='ignore',
                                               duplicate_sequence_action='ignore',
                                               invalid_sequence_action='ignore'))
        seconds = time.perf_counter() - t0
        n_records = len(data)
        n_residues = sum(map(len, data.values()))

    elif mode == 'stream':
        for _header, seq in protfasta.read_fasta_stream(filename):
            n_records += 1
            n_residues += len(seq)

    elif mode == 'stream-unique':
        for _header, seq in protfasta.read_fasta_stream(filename, expect_unique_header=True, silence_warnings=True):
            n_records += 1
            n_residues += len(seq)

    elif mode == 'stream-dedup':
        for _header, seq in protfasta.read_fasta_stream(filename, duplicate_sequence_action='remove', silence_warnings=True):
            n_records += 1
            n_residues += len(seq)

    else:
        raise ValueError('unknown mode %r' % (mode,))

    if not mode.startswith('read_fasta'):
        seconds = time.perf_counter() - t0

    return {
        'mode': mode,
        'seconds': seconds,
        'peak_rss_bytes': _peak_rss_bytes(),
        'n_records': n_records,
        'n_residues': n_residues,
    }


# ---------------------------------------------------------------------------
# driver side
# ---------------------------------------------------------------------------
def _spawn_worker(mode: str, filename: str) -> dict[str, Any]:
    """Run one measurement in a fresh interpreter and return its JSON result."""
    cmd = [sys.executable, os.path.abspath(__file__), '--worker', mode, filename]
    env = dict(os.environ, PYTHONPATH=REPO_ROOT)
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError('worker failed (%s on %s):\n%s' % (mode, filename, proc.stderr))
    return json.loads(proc.stdout.strip().splitlines()[-1])


def _ensure_file(data_dir: str, n_records: int, seed: int) -> str:
    """Return the path of the synthetic file for *n_records*, generating it if needed."""
    from generate_fasta import generate_fasta

    os.makedirs(data_dir, exist_ok=True)
    filename = os.path.join(data_dir, 'synthetic_%d.fasta' % (n_records,))
    if not os.path.exists(filename):
        print('[generate] %s (%d records) ...' % (filename, n_records), end=' ', flush=True, file=sys.stderr)
        t0 = time.perf_counter()
        generate_fasta(filename, n_records, seed=seed)
        print('%.1f s, %.1f MB' % (time.perf_counter() - t0, os.path.getsize(filename) / 1e6), file=sys.stderr)
    return filename


def _warm_cache(filename: str) -> None:
    """Read *filename* sequentially so every mode sees a warm page cache.

    Generating a multi-gigabyte file leaves the OS busy flushing it to disk
    and the page cache in an unknown state, so the first mode timed after
    generation could be several times slower than the rest purely because
    of I/O. A plain sequential read levels the field.
    """
    with open(filename, 'rb') as fh:
        while fh.read(64 << 20):
            pass


def _fmt_count(n: int) -> str:
    return '{:,}'.format(n)


def _fmt_mem(nbytes: float) -> str:
    if nbytes < 0:
        nbytes = 0.0
    if nbytes >= 1e9:
        return '%.2f GB' % (nbytes / 1e9)
    return '%.0f MB' % (nbytes / 1e6)


def _markdown_table(rows: list[dict[str, Any]]) -> str:
    header = ['records', 'file size', 'mode', 'time', 'records/s', 'MB/s', 'peak memory', 'bytes/record']
    lines = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        lines.append('| %s | %s | `%s` | %.2f s | %s | %.0f | %s | %s |' % (
            _fmt_count(r['n_records']),
            _fmt_mem(r['file_bytes']),
            r['mode'],
            r['seconds'],
            _fmt_count(int(r['n_records'] / r['seconds'])) if r['seconds'] > 0 else '-',
            r['file_bytes'] / 1e6 / r['seconds'] if r['seconds'] > 0 else 0.0,
            _fmt_mem(r['memory_bytes']),
            _fmt_count(int(r['memory_bytes'] / r['n_records'])) if r['n_records'] else '-',
        ))
    return '\n'.join(lines)


def _machine_info() -> dict[str, str]:
    import protfasta

    cpu = platform.processor() or platform.machine()
    if sys.platform == 'darwin':
        try:
            cpu = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'], capture_output=True, text=True).stdout.strip() or cpu
        except OSError:
            pass
    return {
        'timestamp': datetime.datetime.now().isoformat(timespec='seconds'),
        'platform': platform.platform(),
        'cpu': cpu,
        'python': platform.python_version(),
        'protfasta': protfasta.__version__,
    }


def run_benchmark(
    sizes: list[int],
    modes: list[str],
    repeats: int,
    data_dir: str,
    results_dir: Optional[str],
    seed: int,
    progress: Callable[[str], None] = lambda msg: print(msg, file=sys.stderr),
) -> list[dict[str, Any]]:
    """Run the full benchmark matrix and return one result row per (size, mode)."""
    info = _machine_info()
    progress('protfasta %s | %s | Python %s | %s' % (info['protfasta'], info['cpu'], info['python'], info['platform']))

    # a dummy file for the import-only baseline (contents irrelevant)
    baseline = _spawn_worker('import-only', __file__)['peak_rss_bytes']
    progress('baseline interpreter + import: %s' % _fmt_mem(baseline))

    rows: list[dict[str, Any]] = []
    for n in sizes:
        filename = _ensure_file(data_dir, n, seed)
        file_bytes = os.path.getsize(filename)
        _warm_cache(filename)
        for mode in modes:
            best: Optional[dict[str, Any]] = None
            for _ in range(repeats):
                res = _spawn_worker(mode, filename)
                if best is None:
                    best = res
                else:
                    best['seconds'] = min(best['seconds'], res['seconds'])
                    best['peak_rss_bytes'] = max(best['peak_rss_bytes'], res['peak_rss_bytes'])
            assert best is not None
            row = {
                'n_records': best['n_records'],
                'n_residues': best['n_residues'],
                'file_bytes': file_bytes,
                'mode': mode,
                'seconds': best['seconds'],
                'peak_rss_bytes': best['peak_rss_bytes'],
                'memory_bytes': max(0, best['peak_rss_bytes'] - baseline),
            }
            rows.append(row)
            progress('  %-20s %-18s %7.2f s  %s' % (_fmt_count(n), mode, row['seconds'], _fmt_mem(row['memory_bytes'])))

    table = _markdown_table(rows)
    print()
    print(table)

    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        stem = os.path.join(results_dir, 'benchmark_%s' % (info['timestamp'].replace(':', '-'),))
        with open(stem + '.json', 'w') as fh:
            json.dump({'machine': info, 'baseline_rss_bytes': baseline, 'repeats': repeats, 'results': rows}, fh, indent=2)
        with open(stem + '.md', 'w') as fh:
            fh.write('# protfasta benchmark\n\n')
            fh.write('- protfasta %s\n- %s\n- Python %s\n- %s\n- %s\n- repeats: %d\n\n' % (
                info['protfasta'], info['cpu'], info['python'], info['platform'], info['timestamp'], repeats))
            fh.write(table + '\n')
        progress('results written to %s.{json,md}' % (stem,))

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sizes', type=int, nargs='+', default=DEFAULT_SIZES,
                        help='record counts to benchmark (default: %s)' % ' '.join(map(str, DEFAULT_SIZES)))
    parser.add_argument('--modes', nargs='+', default=MODES, choices=MODES, help='read modes to benchmark (default: all)')
    parser.add_argument('--repeats', type=int, default=1, help='runs per measurement; fastest time / largest memory kept (default 1)')
    parser.add_argument('--data-dir', default=os.path.join(HERE, 'data'), help='where generated FASTA files are cached')
    parser.add_argument('--results-dir', default=os.path.join(HERE, 'results'), help="where results are written ('' to disable)")
    parser.add_argument('--seed', type=int, default=0, help='seed for the synthetic files (default 0)')
    parser.add_argument('--worker', nargs=2, metavar=('MODE', 'FILE'), help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker:
        mode, filename = args.worker
        print(json.dumps(_run_mode(mode, filename)))
        return

    run_benchmark(args.sizes, args.modes, args.repeats, args.data_dir, args.results_dir or None, args.seed)


if __name__ == '__main__':
    main()
