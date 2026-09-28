# phader

**Peptide insert counting for T7 phage display sequencing data.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

phader turns raw sequencing reads from a phage display biopanning experiment into a ranked
table of peptide counts. It reads FASTQ or FASTA files (optionally gzip-compressed), translates
each read in all six frames, locates the variable insert by the constant leader sequence of the
T7 display construct, applies the filters appropriate to the library design, and tallies
identical peptides per sample.

It is built for the practical question a biopanning experiment asks — *which peptides became
more abundant between rounds of selection* — and for the file sizes that question comes in:
data are read from a directory on the machine running the app, so multi-gigabyte FASTQ files
are never uploaded through a browser.

A graphical interface is provided; the counting engine (`peptide_parser.py`) has no
dependencies beyond the Python standard library and can be imported directly.

## Reproducibility notice

**v1.0.0 is the version used for the analysis reported in the accompanying manuscript.** Its
counting behaviour is frozen deliberately: known rough edges are documented in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md) rather than fixed, so that the results in the paper can be
reproduced exactly from this tag. Behavioural improvements are collected for the next major
version (see [CHANGELOG.md](CHANGELOG.md)).

## Installation

Python 3.9 or newer.

```bash
git clone https://github.com/KarlisPleiko/phader.git
cd phader
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running

```bash
streamlit run streamlit_app.py
```

Streamlit prints a local URL and opens it in your browser. Then:

1. Enter the **absolute path** of the directory holding your FASTQ/FASTA files. The path is
   resolved on the machine running the app, not in the browser.
2. Adjust the **glob patterns** if needed (default
   `*.fastq,*.fq,*.fastq.gz,*.fa,*.fasta`) and press **Scan folder**.
3. Tick the files to analyse — all of them by default.
4. Set the options for your library design (see below).
5. Press **Run parser**, then download the per-sample result tables individually or as a
   single ZIP.

Each input file produces its own result table; nothing is written back to the data directory.

The interface has been exercised end to end on Python 3.13 with Streamlit 1.64.0.

Try it first on the 10-read synthetic dataset in [`examples/`](examples/), which has a known
expected output and takes a second to run.

## Options

Defaults are chosen for a cysteine-constrained library sequenced on an Illumina instrument.

### Library design

| Option | Default | Effect |
|---|---|---|
| Input is FASTA | off | Read FASTA instead of 4-line FASTQ records. Quality filtering is unavailable. |
| Allow stops within first n aa | 0 | Permit stop codons among the first n residues of the insert instead of truncating at the first one. Retained stops appear as `_` in the output. |
| Do NOT require C...C | off | Remove the cysteine requirement **entirely** and report the whole insert up to the first stop. Use this for linear libraries or mixed pools. |
| Require CX4C / CX7C / CX8C / CX10C | off | Require a cysteine, exactly 4/7/8/10 residues, and a second cysteine. These are exact spacings, not minima. |
| Require X7-TEV-His6 | off | Require 7 variable residues followed by a TEV site and a His tag; only the 7 residues are reported. |
| NNK filter | off | Discard reads whose insert codons are incompatible with NNK degeneracy (any third base other than G or T, with the fixed `TGC` cysteine codon exempt). |
| Only forward strand / Only reverse strand | off | Restrict the six-frame search to one strand. Setting both cancels both. |

With no cysteine option selected, the default requirement is a cysteine, at least 3 residues,
and a further cysteine — matching CX4C through CX10C libraries and anything in between.

### Read filtering

| Option | Default | Effect |
|---|---|---|
| Min Phred score (allow up to 3 below) | 0 (off) | Discard a read if **more than three** bases of its insert region fall below this Phred score. On real data this is the single most effective way to suppress error-derived peptide variants. |
| Peptide list to include/exclude | none | A text file of peptide sequences, one per line. By default only these peptides are reported; `X` acts as a wildcard (expanded over all 20 amino acids, up to five per entry). |
| Exclude peptides from peplist | off | Invert the list: report everything except the listed peptides. |

### DNA-level output

| Option | Default | Effect |
|---|---|---|
| Include DNA variants in output | off | Report, under each peptide, every distinct nucleotide sequence that produced it with its own count. The sequence spans the insert through to the end of the read. |
| Spaces between codons in DNA output | off | Split those sequences into space-separated codons. Implies the option above. |
| Etalon peptide list | none | Reference peptide list that enables clustering of DNA variants, merging likely sequencing errors into their parent sequence. |
| Etalon DNA edit distance cutoff | 4.0 | Distance below which a variant is merged into its reference. |
| Don't use dynamic etalon filter | off | Only sequences from the reference list may act as cluster centres. Faster. |
| Remove too-similar DNA | off | Discard merged variants instead of adding their counts to the parent sequence. |

Clustering only engages for variants above a hardcoded 200-copy abundance floor, so it does
nothing on shallowly sequenced libraries. See [KNOWN_ISSUES.md](KNOWN_ISSUES.md) #6.

### Output and performance

| Option | Default | Effect |
|---|---|---|
| Limit output table length | 0 (no limit) | Keep only the n most abundant peptides in the table. Header totals still describe the full result set. |
| CPU workers | 1 | Number of worker processes. Counts are unaffected by this setting; only speed changes. **Raise it for large files** — the default is 1 regardless of available cores. |

## Output format

One whitespace-padded text table per input file, peptides in descending count order:

```
"Unique peptides"  5
 "Total peptides"   7

"Total unique peptides:" 5
"Sum total peptides:" 7
"Peptides"	 "synthetic_demo"
 CSPWNRGVC          3
 CKEFTYHAC          1
```

- **Unique peptides / Total peptides** — the number of distinct peptides and the number of
  reads assigned to them, for this file.
- **Peptide rows** — peptide sequence and its read count. Stop codons, if allowed, appear
  as `_`.
- Fields are padded to fixed widths rather than tab-separated; parse by splitting on
  whitespace.

Counts are deterministic, but peptides **tied on the same count** are not written in a fixed
order, so two runs of the same analysis can produce files that `diff` reports as different
while every count agrees. Sort before comparing, or set `PYTHONHASHSEED=0` — see
[KNOWN_ISSUES.md](KNOWN_ISSUES.md) #17.

With DNA output enabled, each peptide row is followed by its nucleotide variants and counts;
merged variants additionally carry their translation and their distance to the most abundant
variant. A full description of the algorithm, including the filter precedence and the distance
metric, is in [docs/ALGORITHM.md](docs/ALGORITHM.md).

## Performance

Measured on an 8-core laptop (macOS, Python 3.13) using a 223 MB FASTQ file of 468,367
Illumina reads of 201 nt:

| CPU workers | Wall time | Throughput |
|---|---|---|
| 1 | 71 s | ~6,600 reads/s |
| 4 | 20 s | ~23,000 reads/s |

Counts were identical at both settings. Memory use is modest — reads are streamed, and only
the peptide tally is held in memory — but enabling DNA-variant output on a large file stores
every distinct nucleotide sequence and will increase memory use substantially.

## Limitations

The leader sequence that anchors the insert search is specific to the display construct used
here and is hardcoded; other vectors require a source edit. Reads containing ambiguity
characters are counted with the affected codon dropped rather than being rejected, and the
tool does not report why reads were discarded. These and other items are described with their
workarounds in [KNOWN_ISSUES.md](KNOWN_ISSUES.md) — please read it before interpreting counts.

## How to cite

If phader contributes to work you publish, please cite the archived release:

> Pleiko, K. phader: peptide insert counting for T7 phage display sequencing data (v1.0.0).
> Zenodo. DOI: *to be added on release*

Machine-readable metadata is in [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).

## Contact

Questions and bug reports are best raised as
[GitHub issues](https://github.com/KarlisPleiko/phader/issues).
