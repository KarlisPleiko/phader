# Setup checklist

Everything that still needs a human. This file is a scratchpad for the first release — delete
it once the repository is live.

## 1. Fill in the placeholders

| File | Placeholder | What to put there |
|---|---|---|
| `CITATION.cff` | `date-released` | the date you cut the GitHub release |
| `CITATION.cff` | `doi` | the Zenodo DOI (step 4) |
| `README.md` | "How to cite" block | the same DOI, and the manuscript citation once it has one |
| `README.md` | Contact section | the email address you want issues sent to |

Nothing else is templated: owner (`KarlisPleiko`), repository name (`phader`), author
(Karlis Pleiko, ORCID 0000-0002-1073-197X) and the MIT licence are already filled in.

## 2. Check before the first push

- [ ] Confirm the version in `pyproject.toml` (`1.0.0`) and `CITATION.cff` agree with the tag.
- [ ] Confirm the two source files are the exact ones used for the manuscript analysis.
      Their sha256 sums as published here are:
      `peptide_parser.py` `0793d9e66722794f0f961d06106c5d22c3c6421732b7f5de309920fb044811f0`
      `streamlit_app.py` `b54ab636f19bc6c00b64478341139fdfd356f2de5c2a6804d4092b784e6ddad1`
      (`shasum -a 256 peptide_parser.py streamlit_app.py`)
- [ ] Confirm the option settings quoted in the README's worked example match what you
      actually ran for the manuscript. Adjust the example if not.
- [ ] Confirm no sequencing data is staged: `git status --short` should list only text files.
      `.gitignore` already excludes `*.fastq`, `*.fastq.gz`, `data/` and `results/`.

## 3. Create the repository

```bash
cd phader
git init
git add .
git commit -m "phader v1.0.0: peptide insert counting for T7 phage display data"
git branch -M main
git remote add origin git@github.com:KarlisPleiko/phader.git
git push -u origin main
git tag -a v1.0.0 -m "Version used for the manuscript analysis"
git push origin v1.0.0
```

On GitHub, set the repository description and topics (`phage-display`, `t7-phage`,
`biopanning`, `peptide-library`, `fastq`, `ngs`) so the tool is findable.

## 4. Mint a DOI

Journals increasingly require an archived, versioned copy rather than a bare GitHub link.

1. Sign in to [zenodo.org](https://zenodo.org) with your GitHub account.
2. Under *GitHub* in your Zenodo profile, switch the `phader` repository **on**.
3. On GitHub, go to *Releases* → *Draft a new release*, choose the existing `v1.0.0` tag,
   title it `phader v1.0.0`, and publish. Zenodo archives the tag and mints a DOI.
4. Copy the DOI into `CITATION.cff` and the README, commit, and push. (This edit lands after
   the release, which is expected; the DOI resolves to the archived snapshot regardless.)
5. Zenodo also issues a concept DOI that always points at the newest version — cite the
   version DOI in the manuscript so readers land on exactly the code you used.

## 5. Manuscript text

Two places usually need it:

- **Code availability**: name the repository URL, the tag, the DOI and the licence.
- **Methods**: the parameter values you used. `docs/ALGORITHM.md` is written to be
  citable for the algorithm itself so the Methods section can stay short.

## 6. Second app

The differential-homing app is a separate repository. Keep this same structure so the two
look like a family: README with a worked example, `docs/ALGORITHM.md`, `KNOWN_ISSUES.md`,
`CITATION.cff` with the same author metadata, MIT licence, and a Zenodo DOI.
