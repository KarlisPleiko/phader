# Known issues and limitations (v1.0.0)

v1.0.0 is the version used for the analysis in the accompanying manuscript. Its counting
behaviour is therefore **deliberately frozen**: the items below are documented rather than
fixed, so that results reported in the paper can be reproduced exactly from this tag. They
are the work list for the next release.

Each entry states the user-visible effect and, where one exists, the workaround.

## Sequence handling

**1. Codons containing `N` are dropped, not flagged.**
Translation skips any codon that is not in the standard genetic code table, which includes
every codon containing an ambiguity character. The read is still counted, but:

- an `N` inside the insert produces a peptide one residue shorter than the true insert,
  which is tallied as a separate, spurious low-abundance variant;
- an `N` upstream of the insert leaves the peptide correct but shifts the DNA sequence
  reported by *Include DNA variants* by 3 nt per dropped codon.

Both effects were confirmed on synthetic reads. The quality filter does not reliably remove
these reads, because a single `N` is only one sub-threshold base and up to three are
tolerated. If DNA-level output matters, pre-filter reads containing `N` before analysis.

**2. The constant leader is hardcoded.**
The anchor is the 7-residue leader `MLGDPNS` immediately upstream of the variable insert
(exact match, or the 6-residue window `MLGDPN` with up to two substitutions followed by
`S`). Libraries built on a different display vector need a source edit in
`peptide_parser.py` (`mlg_form1` / `mlg_form2`, and the constant `21` = leader length in nt).

**3. `Do NOT require C...C` disables the scaffold check entirely.**
It does not relax the cysteine requirement; no scaffold pattern is applied at all, and the
reported peptide is the whole stop-truncated candidate downstream of the leader. Expect
longer peptides and a larger unique count than with any CXnC option.

**4. The CXnC options are exact spacings, not minima.**
`Require CX7C` matches C + exactly 7 residues + C. A CX7C option run on a CX4C library
returns an empty table rather than a warning.

**5. `Require X7-TEV-His6` requires five histidines.**
The pattern is `[H]{5}` followed by any residues, so a His5 cassette also matches. Only the
seven variable residues are reported, not the TEV site or the tag.

**6. Etalon DNA clustering is inactive at low sequencing depth.**
A DNA variant is promoted to a reference only once its count exceeds a hardcoded 200 copies,
and a candidate is only compared against references more than 200x more abundant than
itself. Below those thresholds the option runs and silently changes nothing.

## Reporting

**7. Per-read rejection statistics are computed and then discarded.**
The engine tracks how far each read progressed (no leader found, leader but no scaffold
match, quality rejection, and so on) but nothing surfaces those counters. There is currently
no way to see *why* reads were dropped; only accepted reads appear in the totals.

**8. One count column per output file.**
The table format supports several columns, but each input file is written to its own results
file, so cross-round comparisons must be assembled outside the tool (join the per-file
tables on the peptide column).

**9. Output is whitespace-padded text, not strict TSV.**
Fields are space-padded to fixed widths and the header labels are quoted. Parse with
`str.split()` on whitespace rather than splitting on tabs.

**10. `Limit output table length` truncates rows only.**
The unique- and total-peptide figures in the header always describe the full result set, not
the truncated table.

## Reproducibility

**17. The order of peptides with identical counts varies between runs.**
Counts themselves are fully deterministic — repeated runs of the same file with the same
options always produce the same numbers, at any worker setting. Rows are sorted by count
only, and peptides tied on count are emitted in Python set-iteration order, which depends on
the per-process string hash seed. Two runs of the same analysis therefore produce files that
`diff` reports as different even though every count agrees.

This matters when comparing result files byte-for-byte or checksumming them for an archive.
Either sort the table before comparing (`sort -k2,2nr -k1,1`), or set the hash seed:

```bash
PYTHONHASHSEED=0 streamlit run streamlit_app.py
```

`examples/expected_output_default.txt` was generated with `PYTHONHASHSEED=0` and is
reproducible byte-for-byte under that setting.

## Interface and performance

**11. The CPU-worker slider defaults to 1.**
The default is pinned to one worker regardless of how many cores the machine has. Move the
slider up manually for large files; counts are identical at any worker setting (verified).

**12. The progress bar does not track progress.**
It jumps from 0 to 100 % when the whole batch finishes. Long runs look stalled; they are not.

**13. The engine keeps configuration in module-level state.**
Two analyses running concurrently in the same Streamlit server process can interfere. Run one
analysis at a time, and treat the app as single-user.

**14. macOS uses the `fork` start method.**
It is set explicitly at import so that workers appear as real processes under Streamlit.
This is the start method Apple documents as unsafe in general; it is safe here because
workers touch only pure-Python state, but it is a reason to keep the worker count modest.

**15. Deprecated Streamlit argument.**
`use_container_width` is deprecated in recent Streamlit releases. `requirements.txt` therefore
pins `streamlit<2`.

**16. The `report` option is accepted and ignored.**
It is retained only so the function signature matches the original command-line script.
