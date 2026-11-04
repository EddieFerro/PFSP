# PFSP: Prefix-Free Syncmer Parsing

PFSP parses a collection of sequences (for example a pangenome of haplotype assemblies) into a dictionary of
distinct phrases and a parse, like prefix-free parsing (PFP), but with phrase boundaries set by closed syncmers
instead of a rolling-hash window. This bounds the size of every phrase.

## Build

Requirements: a C++17 compiler, CMake ≥ 3.20, zlib. Other dependencies (CLI11, kseq, MurmurHash3) are downloaded
at configure time at pinned versions.

```bash
git clone https://github.com/<user>/PFSP.git
cd PFSP
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j        # -> build/PFSP
```

## Usage

```bash
build/PFSP input.fa.gz -o out/prefix -k 200 -s 20 --seed 0
```

| Option | Meaning |
|---|---|
| `input` | FASTA/FASTQ, optionally gzipped (`/dev/stdin` to stream) |
| `-o, --output` | output prefix (default: the input path) |
| `-k, --syncmer-length` | syncmer length K (default 200) |
| `-s, --smer-length` | s-mer length S, at most 21 (default 20) |
| `--seed` | seed for hashing (default 0) |

Outputs:

| File | Contents |
|---|---|
| `<prefix>.dict` | the distinct syncmers (phrases), sorted, 3 bits per base |
| `<prefix>.counts` | occurrences of each phrase (uint32, same order) |
| `<prefix>.parse` | the parse: 1-based phrase ranks (uint32) |
| `<prefix>.gaps` | distance from the previous phrase start (uint16), one per parse entry |
| `<prefix>.meta` | parameters and, per record, its length and range in the parse |

Each record is padded with S `$` symbols at both ends, so it starts and ends with a phrase; records shorter than
K − 2S are skipped (listed in `.meta` with 0 entries). Consecutive phrases overlap, and the gap between them is at
most K − S + 1, which bounds the phrase size.
The same `-k`, `-s` and `--seed` must be used wherever the parse is used later (for example by SIFT).

## Citation

<!-- paper reference -->

## License

GNU General Public License v3.0 (see `LICENSE`).
