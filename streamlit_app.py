# streamlit_app.py
import os
import glob
import pathlib
import tempfile
from zipfile import ZipFile

import streamlit as st
from peptide_parser import run_on_files

st.set_page_config(page_title="FASTQ → Peptide Counter (Local Folder Mode)", layout="wide")

st.title("FASTQ → Peptide Counter")
st.caption("Point to a local folder (server-side) with FASTQ/FASTA files (incl. .gz). No browser upload of big files.")

# ---------------- Options ----------------
with st.expander("Options", expanded=True):
    col1, col2, col3 = st.columns(3)
    with col1:
        fasta = st.checkbox("Input is FASTA (default FASTQ)", value=False)
        stops = st.number_input("Allow stops within first n aa", min_value=0, max_value=50, value=0, step=1)
        nocys = st.checkbox("Do NOT require C...C", value=False)
        his = st.checkbox("Require X7-TEV-His6", value=False)
        qual = st.number_input("Min Phred score (allow up to 3 below)", min_value=0, max_value=60, value=0, step=1)
        nnk = st.checkbox("NNK filter", value=False)
        limit = st.number_input("Limit output table length (0 = no limit)", min_value=0, value=0, step=100)
    with col2:
        cx4c = st.checkbox("Require CX4C", value=False)
        cx7c = st.checkbox("Require CX7C", value=False)
        cx8c = st.checkbox("Require CX8C", value=False)
        cx10c = st.checkbox("Require CX10C", value=False)
        rev = st.checkbox("Only reverse strand", value=False)
        fwd = st.checkbox("Only forward strand", value=False)
        dnaspace = st.checkbox("Spaces between codons in DNA output", value=False)
    with col3:
        dna = st.checkbox("Include DNA variants in output", value=False)
        etalon_file = st.file_uploader("Etalon peptide list (optional)", type=["txt","csv","tsv"], accept_multiple_files=False)
        etaloncutoff = st.number_input("Etalon DNA edit distance cutoff", min_value=0.0, max_value=20.0, value=4.0, step=0.5)
        oldetalon = st.checkbox("Don't use dynamic etalon filter (faster)", value=False)
        etalonremove = st.checkbox("Remove too-similar DNA (instead of marking)", value=False)
        mp_workers = st.slider("CPU workers", 1, os.cpu_count() or 4, min(1, os.cpu_count() or 4))

    peplist_file = st.file_uploader("Peptide list to include/exclude (optional)", type=["txt","csv","tsv"], accept_multiple_files=False)
    excludepl = st.checkbox("Exclude peptides from peplist (otherwise include only those)", value=False)

st.markdown("---")

# ---------------- Local folder mode ----------------
st.subheader("Local folder (no upload)")

# Persist scanned files across reruns
if "scanned_files" not in st.session_state:
    st.session_state.scanned_files = []

root = st.text_input(
    "Absolute path to a local folder on the server",
    value=str(pathlib.Path.cwd())
)
patterns = st.text_input(
    "Glob patterns (comma-separated)",
    value="*.fastq,*.fq,*.fastq.gz,*.fa,*.fasta"
)

def scan_folder(folder, pats):
    files = []
    for p in pats:
        files.extend(glob.glob(os.path.join(folder, p)))
    files = sorted(set(files))
    return files

if st.button("Scan folder"):
    pats = [p.strip() for p in patterns.split(",") if p.strip()]
    files = scan_folder(root, pats)
    st.session_state.scanned_files = files
    if not files:
        st.warning("No files found for given patterns.")
    else:
        st.success(f"Found {len(files)} files.")

files = st.session_state.scanned_files

if files:
    st.write("Select files to process:")
    # “Select all” convenience
    select_all = st.checkbox("Select all", value=True)
    default_selection = files if select_all else []
    selected = st.multiselect("Files", files, default=default_selection)
else:
    selected = []

run_btn = st.button("Run parser", type="primary", use_container_width=True)

# ---------------- Run ----------------
if run_btn:
    if not selected:
        st.error("No files selected.")
        st.stop()

    # Save small aux files (etalon, peplist) to temp paths if provided
    tmpdir = tempfile.mkdtemp(prefix="pepc_")

    etalon_path = None
    if etalon_file is not None:
        etalon_path = os.path.join(tmpdir, "etalon.txt")
        with open(etalon_path, "wb") as f:
            f.write(etalon_file.read())

    peplist_path = None
    if peplist_file is not None:
        peplist_path = os.path.join(tmpdir, "peplist.txt")
        with open(peplist_path, "wb") as f:
            f.write(peplist_file.read())

    st.info(f"Processing {len(selected)} file(s) from disk…")
    progress = st.progress(0)

    results_map = run_on_files(
        selected,
        fasta=fasta,
        stops=stops,
        nocys=nocys,
        rev=rev,
        fwd=fwd,
        limit=limit,
        peplist_path=peplist_path,
        excludepl=excludepl,
        dna=dna,
        dnaspace=dnaspace,
        etalon_path=etalon_path,
        etaloncutoff=etaloncutoff,
        oldetalon=oldetalon,
        etalonremove=etalonremove,
        nnk=nnk,
        his=his,
        qual=qual,
        cx7c=cx7c,
        cx8c=cx8c,
        cx10c=cx10c,
        cx4c=cx4c,
        report=False,
        mp_workers=mp_workers,
    )
    progress.progress(100)

    # Per-file downloads + a ZIP bundle
    zip_path = os.path.join(tmpdir, "all_results.zip")
    with ZipFile(zip_path, "w") as zf:
        for base, text in results_map.items():
            resname = f"{base}_results.txt"
            respath = os.path.join(tmpdir, resname)
            with open(respath, "w") as f:
                f.write(text)
            zf.write(respath, arcname=resname)

    st.success("Done. Download results below:")
    for base, text in results_map.items():
        fname = f"{base}_results.txt"
        st.download_button(
            label=f"Download {fname}",
            data=text.encode("utf-8"),
            file_name=fname,
            mime="text/plain",
        )

    with open(zip_path, "rb") as f:
        st.download_button(
            label="Download ALL results as ZIP",
            data=f.read(),
            file_name="all_results.zip",
            mime="application/zip",
        )

st.markdown("---")
st.caption("Tip: works best when the data directory is on a fast local disk or mounted volume. Supports .fastq, .fq, .fa, .fasta, and .fastq.gz.")
