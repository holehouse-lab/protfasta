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
| `stream` | `protfasta.read_fasta_stream(f)` with default options, fully consumed. This is the flat-memory configuration. |
| `stream-unique` | `read_fasta_stream(f, expect_unique_header=True)` - keeps a running set of every header. |
| `stream-dedup` | `read_fasta_stream(f, duplicate_sequence_action='remove')` - keeps a running set of 16-byte sequence digests. |

Every (size, mode) measurement runs in a fresh subprocess so peak memory is not contaminated by earlier runs. **Time** covers only the read call itself (not interpreter start-up or imports). **Memory** is the subprocess's peak resident set size (`ru_maxrss`) minus the peak of a subprocess that only imports protfasta (about 24 MB), i.e. the memory the parse itself needed. For `read_fasta` that includes the returned data structure, which is the whole point of the comparison with the streaming modes. With `--repeats N` the fastest time and the largest memory over N runs are kept.

## Reference results

Apple M3 Max, macOS 26.5, Python 3.12.13, protfasta 0.1.24 (August 2026). `read_fasta` bytes/record is the retained size of the result: roughly the sequence itself (~420 residues) plus ~400 bytes of Python object overhead (header string, sequence string, dictionary entry).

- protfasta 0.1.23+2.gbdc9d94.dirty
- Apple M3 Max
- Python 3.12.13
- macOS-26.5.1-arm64-arm-64bit
- 2026-08-29T23:13:27
- repeats: 2

| records | file size | mode | time | records/s | MB/s | peak memory | bytes/record |
|---|---|---|---|---|---|---|---|
| 10,000 | 5 MB | `read_fasta` | 0.02 s | 489,307 | 255 | 8 MB | 822 |
| 10,000 | 5 MB | `read_fasta-nocheck` | 0.01 s | 691,479 | 361 | 8 MB | 794 |
| 10,000 | 5 MB | `stream` | 0.02 s | 557,439 | 291 | 0 MB | 39 |
| 10,000 | 5 MB | `stream-unique` | 0.02 s | 513,700 | 268 | 2 MB | 208 |
| 10,000 | 5 MB | `stream-dedup` | 0.03 s | 375,117 | 196 | 3 MB | 278 |
| 100,000 | 52 MB | `read_fasta` | 0.22 s | 451,177 | 236 | 84 MB | 838 |
| 100,000 | 52 MB | `read_fasta-nocheck` | 0.16 s | 628,699 | 329 | 81 MB | 813 |
| 100,000 | 52 MB | `stream` | 0.18 s | 545,727 | 286 | 0 MB | 2 |
| 100,000 | 52 MB | `stream-unique` | 0.20 s | 490,373 | 257 | 22 MB | 216 |
| 100,000 | 52 MB | `stream-dedup` | 0.28 s | 362,627 | 190 | 31 MB | 314 |
| 1,000,000 | 527 MB | `read_fasta` | 2.35 s | 425,359 | 224 | 813 MB | 812 |
| 1,000,000 | 527 MB | `read_fasta-nocheck` | 1.72 s | 582,242 | 307 | 799 MB | 798 |
| 1,000,000 | 527 MB | `stream` | 1.81 s | 553,206 | 292 | 0 MB | 0 |
| 1,000,000 | 527 MB | `stream-unique` | 2.06 s | 485,016 | 256 | 211 MB | 210 |
| 1,000,000 | 527 MB | `stream-dedup` | 2.96 s | 338,262 | 178 | 293 MB | 293 |
| 10,000,000 | 5.30 GB | `read_fasta` | 24.75 s | 404,058 | 214 | 7.93 GB | 793 |
| 10,000,000 | 5.30 GB | `read_fasta-nocheck` | 18.93 s | 528,310 | 280 | 7.84 GB | 783 |
| 10,000,000 | 5.30 GB | `stream` | 19.96 s | 500,948 | 265 | 0 MB | 0 |
| 10,000,000 | 5.30 GB | `stream-unique` | 24.16 s | 413,843 | 219 | 1.99 GB | 198 |
| 10,000,000 | 5.30 GB | `stream-dedup` | 28.62 s | 349,434 | 185 | 2.76 GB | 276 |

Takeaways:

* **Parsing runs at roughly 200-300 MB/s, or 400,000-550,000 records/s**, in pure Python, and scales linearly. A 100,000-record proteome (~50 MB) loads in a fifth of a second; a million records in ~2 seconds; ten million (5.3 GB) in 20 seconds streamed or 25 seconds fully loaded.
* **`read_fasta` memory is ~800 bytes per record**, i.e. about 1.5x the size of the file on disk. Ten million records need ~8 GB. If that is a problem, `read_fasta_stream` does the same parse in constant memory: it never rose above the measurement floor at any size.
* **The opt-in streaming checks cost ~200 bytes/record** (`expect_unique_header=True`, a set of header strings) **to ~280 bytes/record** (`duplicate_sequence_action='remove'`, a set of 16-byte digests plus set overhead). That is still 3-4x lighter than a full load, but it does grow with the file - which is why they are off by default when streaming.
* The default `read_fasta` checks (header uniqueness and invalid-residue validation) are cheap relative to parsing: about 30% on top of the no-check parse.
* Both `read_fasta` and `read_fasta_stream` start paying off at exactly the file sizes you would hope: there is no per-file overhead worth mentioning, so 10,000 records take 20 ms.

## Notes and caveats

* Synthetic residues are uniformly random, so the sequences compress worse and validate identically to real ones; real headers are often shorter than the UniProt-style ones used here, which slightly reduces the per-record memory in practice.
* Results are single-threaded wall-clock on an otherwise idle machine; disk read time is included but the files are small enough to be in the page cache after generation.
* `ru_maxrss` is reported in bytes on macOS and kilobytes on Linux; the driver normalises this.
