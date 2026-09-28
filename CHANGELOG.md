# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - unreleased

First public release. **This is the version used for the peptide-counting analysis in the
accompanying manuscript**, and its counting behaviour is frozen for reproducibility: the
engine (`peptide_parser.py`) and interface (`streamlit_app.py`) are published exactly as
they were run on the manuscript data.

### Added

- Streamlit interface (`streamlit_app.py`) that analyses a server-side directory of
  FASTQ/FASTA files, with no browser upload of sequencing data.
- Counting engine (`peptide_parser.py`): six-frame translation, leader-anchored insert
  extraction with a mismatch-tolerant fallback, scaffold filters (generic C...C, CX4C, CX7C,
  CX8C, CX10C, X7-TEV-His6), Phred quality filter, NNK codon filter, peptide allow/deny
  lists with `X` wildcard expansion, optional DNA-variant reporting, and edit-distance
  clustering of DNA variants to collapse sequencing errors.
- Support for plain and gzip-compressed inputs and for multiprocessing across CPU cores.
- Documentation: `README.md`, `docs/ALGORITHM.md`, `KNOWN_ISSUES.md`, `CITATION.cff`.

### Known issues

See [KNOWN_ISSUES.md](KNOWN_ISSUES.md). Nothing listed there is fixed in this release by
design; fixes are deferred to keep the manuscript's results reproducible.

## Planned

### [2.0.0]

Refactor for the software-focused publication: configuration passed explicitly instead of
module-level state, per-read rejection statistics surfaced in the output, a combined
multi-file count table, a command-line entry point, a test suite, and corrected handling of
reads containing ambiguity characters. Changes that alter counts will be documented against
the v1.0.0 behaviour.
