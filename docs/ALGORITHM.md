# How phader counts peptides

This document describes exactly what the counting engine does, in the order it does it. It is
the reference for the Methods section of work that uses phader v1.0.0. Behaviour that is
surprising or known to be imperfect is collected separately in
[KNOWN_ISSUES.md](../KNOWN_ISSUES.md).

## Input

One or more FASTQ or FASTA files, plain or gzip-compressed. FASTQ records are read as strict
four-line blocks; quality strings are retained because they are needed by the quality filter.
Reads are processed independently, so they can be distributed across CPU workers without
affecting the result.

## Per-read processing

### 1. Six-frame translation

Each read is examined in six frames — reading-frame offsets 2, 1 and 0 on the reverse
complement, then offsets 2, 1 and 0 on the forward strand. Each frame is translated with the
standard genetic code; stop codons become `_`. *Only reverse strand* and *Only forward strand*
restrict the search to three frames (setting both cancels both).

Frames are searched twice: first requiring an exact leader match, then allowing leader
mismatches (step 2). The first frame that yields an acceptable peptide is used and the search
stops, so an exact leader match is always preferred over a mismatched one.

### 2. Anchoring on the constant leader

The variable insert is located by the constant 7-residue leader `MLGDPNS` that immediately
precedes it in the display construct. The first pass requires this sequence exactly. The
second pass accepts the 6-residue window `MLGDPN` with up to two substitutions, followed by
`S`, which recovers inserts from reads carrying sequencing errors in the leader itself.
Everything downstream of the leader becomes the candidate insert.

Reads with no leader match in any frame are not counted.

### 3. Truncation at stop codons

By default the candidate ends at the first stop codon. *Allow stops within first n aa* permits
stop symbols among the first n residues; when this is used, retained stops appear as `_` in
the reported peptide sequence.

### 4. Scaffold requirement

One structural requirement is applied to the candidate, in this precedence order:

| Option | Requirement | Reported peptide |
|---|---|---|
| Require X7-TEV-His6 | 7 residues, then `ENLYFQR` with up to 2 substitutions, then 5+ histidines | the 7 variable residues only |
| Require CX4C | `C` + exactly 4 residues + `C` | from the first `C` to the end of the candidate |
| Require CX7C | `C` + exactly 7 residues + `C` | as above |
| Require CX8C | `C` + exactly 8 residues + `C` | as above |
| Require CX10C | `C` + exactly 10 residues + `C` | as above |
| *(default)* | `C` + at least 3 residues + `C` | from the first `C` to the **last** `C` in the candidate, plus any trailing residues |
| Do NOT require C...C | none | the whole stop-truncated candidate |

Candidates that fail the requirement are not counted.

### 5. Quality filter

Active when *Min Phred score* is greater than 0. Phred scores are decoded from the quality
string (offset 33) over the DNA window corresponding to the accepted insert and downstream
read. The read is discarded if **more than three** bases in that window fall below the
threshold; up to three are tolerated. On reverse-strand hits the quality string is reversed
alongside the sequence, so the window corresponds to the same bases.

### 6. NNK filter

Active when *NNK filter* is selected. Libraries built with NNK degenerate codons carry only
`G` or `T` at every third position. Each codon of the insert is checked and the read is
discarded if any third position is `A` or `C`. The codon `TGC` is exempt, since the flanking
cysteines of the scaffolds are encoded by a fixed `TGC` codon rather than by the degenerate
region.

## Aggregation

Accepted peptides are tallied per input file. When *Include DNA variants in output* is
selected, the engine also records, for each peptide, every distinct nucleotide sequence that
produced it together with its count. The recorded DNA runs from the start of the insert to the
end of the read, so it includes constant 3' vector sequence.

### Optional clustering of DNA variants

Sequencing errors inflate the number of distinct DNA variants observed for one peptide.
Supplying an *etalon peptide list* enables a clustering pass that folds error-derived variants
into their parent sequence:

1. Variants of the listed peptides that exceed 200 copies become reference sequences.
2. Every other variant, processed from most to least abundant, is compared against all
   reference sequences that are more than 200x more abundant than itself.
3. Distance is a modified edit distance: substitutions cost 1.9, insertions and deletions
   cost 1, and runs of consecutive matching positions reduce the score. The weighting
   therefore prefers to explain a difference between two variants as an indel rather than
   as a substitution. This is an unusual choice for Illumina data, where substitution
   errors dominate, and it is reported here as the behaviour of the code rather than as a
   recommendation.
4. A variant within *Etalon DNA edit distance cutoff* of its closest reference is merged into
   it; otherwise it becomes a reference itself once it exceeds 200 copies.

*Don't use dynamic etalon filter* disables step 4's promotion of new references, which is
faster. *Remove too-similar DNA* discards merged variants instead of adding their counts to
the parent.

Because of the hardcoded 200-copy thresholds, this pass has no effect on shallowly sequenced
libraries.

### Peptide allow/deny list

A supplied peptide list restricts the output to the listed peptides, or, with *Exclude
peptides from peplist*, removes them. List entries may contain `X` as a wildcard, expanded
combinatorially over the 20 amino acids (entries with more than five `X` are ignored).

## Output

One whitespace-padded text table per input file, sorted by descending count:

```
"Unique peptides"  5
 "Total peptides"   7

"Total unique peptides:" 5
"Sum total peptides:" 7
"Peptides"	 "sample_name"
 CSPWNRGVC          3
 CKEFTYHAC          1
```

The header reports the number of distinct peptides and the total number of counted reads for
the file; these always describe the complete result set, even when *Limit output table length*
truncates the rows below. With DNA output enabled, each peptide row is followed by its DNA
variants and their counts; for merged variants the row additionally carries the translated
sequence and its distance to the most abundant variant of that peptide.
