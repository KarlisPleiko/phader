# Synthetic example

`synthetic_demo.fastq` contains 10 hand-built reads (no experimental data) that exercise every
branch of the counting engine. It exists so that a new user can confirm an installation
behaves correctly in a few seconds.

Each read follows the construct
`[20 nt 5' flank][leader MLGDPNS][insert][TAA][48 nt 3' flank]`, is 119 nt long, and carries
Phred 40 (`I`) qualities unless stated otherwise.

| # | read id | what it tests | counted by default? |
|---|---|---|---|
| 1 | `fwd_CX7C_CSPWNRGVC` | forward strand, exact leader, insert C-SPWNRGV-C | yes |
| 2 | `rev_CX7C_CKEFTYHAC` | reverse strand, insert C-KEFTYHA-C | yes |
| 3 | `fwd_CX7C_CSPWNRGVC_dup` | duplicate of read 1 | yes |
| 4 | `fwd_internal_stop` | stop codon at insert position 4 | no |
| 5 | `fwd_leader_2mm` | two substitutions in the leader | yes, via the mismatch-tolerant anchor |
| 6 | `fwd_with_N` | one `N` inside the insert | yes, as a truncated variant (see KNOWN_ISSUES #1) |
| 7 | `fwd_lowqual` | six bases at Phred 2 in the insert | yes by default, rejected at Phred >= 20 |
| 8 | `fwd_nonNNK` | insert codon `CCA` (third base not G/T) | yes by default, rejected by the NNK filter |
| 9 | `fwd_CX2C` | only two residues between the cysteines | no |
| 10 | `junk_no_leader` | `(ACGT)x30`, no leader | no |

## Expected result with default options

Run the app, point it at this directory, leave every option at its default and press
**Run parser**. The downloaded `synthetic_demo_results.txt` should carry these counts:

```
"Unique peptides"  5
 "Total peptides"   7

"Total unique peptides:" 5
"Sum total peptides:" 7
"Peptides"	 "synthetic_demo"
 CSPWNRGVC          3
 CKEFTYHAC          1
 CQQWNRGVC          1
 CSPNRGVC           1
 CPPWNRGVC          1
```

The four peptides tied at one count may appear in a different order than shown — see
[KNOWN_ISSUES.md](../KNOWN_ISSUES.md) #17. `expected_output_default.txt` in this directory was
produced with `PYTHONHASHSEED=0` and is reproduced byte-for-byte by

```bash
PYTHONHASHSEED=0 streamlit run streamlit_app.py
```

Reading the table: `CSPWNRGVC` is counted 3 times (reads 1, 3 and 5 — the mismatched leader is
recovered and merged with the exact matches). `CSPNRGVC` is read 6, one residue short because
the codon containing `N` was dropped in translation. `CQQWNRGVC` is read 7 and `CPPWNRGVC` is
read 8; both are admitted by default and removed by the quality and NNK filters respectively.

Useful checks from the same file:

| Change | Expected effect |
|---|---|
| Min Phred score = 20 | `CQQWNRGVC` disappears |
| NNK filter on | `CPPWNRGVC` disappears |
| Do NOT require C...C | `CSPC` and `CSP` appear (7 unique peptides) |
| Allow stops within first 4 aa | `CSP_NRGVC` appears |
| Require CX4C | empty table — every insert here is CX7C |
| Only reverse strand | only `CKEFTYHAC` remains |
| CPU workers = 4 | counts unchanged |
