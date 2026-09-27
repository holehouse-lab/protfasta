# protfasta benchmarks

This directory holds a small, self-contained harness for measuring how protfasta's parsing speed and memory footprint scale with file size, from tens of thousands of records up to the tens of millions.

Two scripts:

* `generate_fasta.py` - writes a synthetic protein FASTA file with a given number of records. Sequence lengths follow a log-normal distribution (median ~350 residues, clipped to 30-5000, so a mean of ~420), residues are uniform over the 20 standard amino acids, headers mimic UniProt and sequences are wrapped at 60 residues per line. Generation is vectorised with numpy, so a 10-million-record (~5.3 GB) file takes about a minute. The same seed always gives the same file.
* `benchmark_protfasta.py` - the driver. For each requested size it generates (and caches) the input file and then times a set of read modes.

## Running it

From the repository root:

```bash
# default: 10k, 100k, 1M and 10M records, all modes, one run each
python benchmark/benchmark_protfasta.py

# smaller / quicker
python benchmark/benchmark_protfasta.py --sizes 10000 100000 1000000 --repeats 3

# just the two headline modes
python benchmark/benchmark_protfasta.py --modes read_fasta stream
```

Generated inputs are cached in `benchmark/data/` and results are written (as JSON and Markdown) to `benchmark/results/`; both directories are git-ignored. Note that if your checkout lives somewhere that is synced (Dropbox, iCloud, ...) you probably want `--data-dir /tmp/protfasta_bench` or similar, since the full set of inputs is ~6 GB. The 10M-record input is ~5.3 GB on disk, and `read_fasta` on it needs ~8 GB of RAM, so drop that size (`--sizes 10000 100000 1000000`) on a small machine. numpy is required for the generator; nothing else beyond protfasta itself.

## What is measured

Modes:

| mode | what it does |
|---|---|
| `read_fasta` | `protfasta.read_fasta(f)` with default options: unique headers enforced, duplicate records fail, invalid residues fail. Returns a dictionary. |
| `read_fasta-nocheck` | `read_fasta` with every check disabled - the raw cost of parsing plus building the result. |
| `read_fasta-dedup` | `read_fasta(f, duplicate_sequence_action='remove')` - the default checks plus duplicate-sequence removal. |
| `stream` | `protfasta.read_fasta_stream(f)` with default options, fully consumed. This is the flat-memory configuration. |
| `stream-unique` | `read_fasta_stream(f, expect_unique_header=True)` - keeps a running set of every header. |
| `stream-dedup` | `read_fasta_stream(f, duplicate_sequence_action='remove')` - keeps a running set of 16-byte sequence digests. |

Every (size, mode) measurement runs in a fresh subprocess so peak memory is not contaminated by earlier runs. **Time** covers only the read call itself (not interpreter start-up or imports, and for the `read_fasta` modes not the pass that counts the residues in the result afterwards). **Memory** is the subprocess's peak resident set size (`ru_maxrss`) minus the peak of a subprocess that only imports protfasta (about 23 MB), i.e. the memory the parse itself needed. For `read_fasta` that includes the returned data structure, which is the whole point of the comparison with the streaming modes. With `--repeats N` the fastest time and the largest memory over N runs are kept.

## Reference results

Apple M3 Max, macOS 26.5, Python 3.12.13, protfasta 0.1.24 (September 2026). `read_fasta` bytes/record is the retained size of the result: roughly the sequence itself (~420 residues) plus ~400 bytes of Python object overhead (header string, sequence string, dictionary entry).

- protfasta 0.1.23+4.g0ccb7f5.dirty
- Apple M3 Max
- Python 3.12.13
- macOS-26.5.1-arm64-arm-64bit
- 2026-09-26T21:36:40
- repeats: 2

| records | file size | mode | time | records/s | MB/s | peak memory | bytes/record |
|---|---|---|---|---|---|---|---|
| 10,000 | 5 MB | `read_fasta` | 0.02 s | 541,742 | 282 | 8 MB | 825 |
| 10,000 | 5 MB | `read_fasta-nocheck` | 0.01 s | 680,876 | 355 | 8 MB | 753 |
| 10,000 | 5 MB | `read_fasta-dedup` | 0.02 s | 511,471 | 267 | 8 MB | 802 |
| 10,000 | 5 MB | `stream` | 0.02 s | 558,168 | 291 | 0 MB | 0 |
| 10,000 | 5 MB | `stream-unique` | 0.02 s | 519,988 | 271 | 2 MB | 190 |
| 10,000 | 5 MB | `stream-dedup` | 0.03 s | 380,745 | 199 | 1 MB | 111 |
| 100,000 | 52 MB | `read_fasta` | 0.19 s | 516,710 | 270 | 82 MB | 817 |
| 100,000 | 52 MB | `read_fasta-nocheck` | 0.15 s | 658,979 | 345 | 81 MB | 807 |
| 100,000 | 52 MB | `read_fasta-dedup` | 0.22 s | 460,688 | 241 | 83 MB | 826 |
| 100,000 | 52 MB | `stream` | 0.18 s | 563,932 | 295 | 0 MB | 0 |
| 100,000 | 52 MB | `stream-unique` | 0.20 s | 507,419 | 266 | 21 MB | 210 |
| 100,000 | 52 MB | `stream-dedup` | 0.27 s | 375,596 | 197 | 13 MB | 132 |
| 1,000,000 | 527 MB | `read_fasta` | 1.97 s | 506,356 | 267 | 814 MB | 813 |
| 1,000,000 | 527 MB | `read_fasta-nocheck` | 1.56 s | 640,946 | 338 | 796 MB | 795 |
| 1,000,000 | 527 MB | `read_fasta-dedup` | 2.22 s | 449,694 | 237 | 831 MB | 831 |
| 1,000,000 | 527 MB | `stream` | 1.75 s | 570,925 | 301 | 0 MB | 0 |
| 1,000,000 | 527 MB | `stream-unique` | 2.01 s | 496,503 | 261 | 211 MB | 210 |
| 1,000,000 | 527 MB | `stream-dedup` | 2.71 s | 368,673 | 194 | 130 MB | 129 |
| 10,000,000 | 5.30 GB | `read_fasta` | 20.90 s | 478,397 | 253 | 7.93 GB | 793 |
| 10,000,000 | 5.30 GB | `read_fasta-nocheck` | 17.00 s | 588,118 | 312 | 7.84 GB | 783 |
| 10,000,000 | 5.30 GB | `read_fasta-dedup` | 23.19 s | 431,185 | 228 | 8.03 GB | 802 |
| 10,000,000 | 5.30 GB | `stream` | 17.47 s | 572,473 | 303 | 0 MB | 0 |
| 10,000,000 | 5.30 GB | `stream-unique` | 19.84 s | 504,063 | 267 | 1.99 GB | 198 |
| 10,000,000 | 5.30 GB | `stream-dedup` | 26.61 s | 375,833 | 199 | 1.18 GB | 117 |

Takeaways:

* **Parsing runs at roughly 250-300 MB/s, or 450,000-600,000 records/s**, in pure Python, and scales linearly. A 100,000-record proteome (~50 MB) loads in under a fifth of a second; a million records in ~2 seconds; ten million (5.3 GB) in 17-18 seconds streamed or 21 seconds fully loaded.
* **`read_fasta` memory is ~800 bytes per record**, i.e. about 1.5x the size of the file on disk. Ten million records need ~8 GB. If that is a problem, `read_fasta_stream` does the same parse in constant memory: it never rose above the measurement floor at any size.
* **The opt-in streaming checks cost ~120 bytes/record** (`duplicate_sequence_action='remove'`, a set of 16-byte digests) **to ~200 bytes/record** (`expect_unique_header=True`, which has to keep every header string). That is still 4-7x lighter than a full load, but it does grow with the file - which is why they are off by default when streaming.
* The default `read_fasta` checks (header uniqueness and invalid-residue validation) are cheap relative to parsing: about 25% on top of the no-check parse. Removing duplicate sequences as well costs another ~10% and almost no memory, because `read_fasta` compares the sequence strings it already holds; the streaming equivalent (`stream-dedup`) has to hash every sequence, since it does not keep them, and pays about 50%.
* Both `read_fasta` and `read_fasta_stream` start paying off at exactly the file sizes you would hope: there is no per-file overhead worth mentioning, so 10,000 records take 20 ms.

## Notes and caveats

* Synthetic residues are uniformly random, so the sequences compress worse and validate identically to real ones; real headers are often shorter than the UniProt-style ones used here, which slightly reduces the per-record memory in practice.
* Results are single-threaded wall-clock on an otherwise idle machine; disk read time is included but the files are small enough to be in the page cache after generation.
* `ru_maxrss` is reported in bytes on macOS and kilobytes on Linux; the driver normalises this.
