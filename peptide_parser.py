# peptide_parser.py
# Multiprocessing-safe + gzip + Streamlit hot-reload safe.

import re
import gzip
from collections import defaultdict
from operator import itemgetter
from os import path
from io import StringIO
import multiprocessing as mp
import types
import platform

# macOS: prefer real processes (fork) so Streamlit shows multiple workers
if platform.system() == "Darwin":
    try:
        mp.set_start_method("fork", force=True)
    except RuntimeError:
        # Already set in this interpreter session
        pass

# ----------------- Globals -----------------
args2 = None            # parent-process config (arg-like)
titles = None
pool = None

# Per-process storage for worker-local args/regex (populated by _init_worker)
_WORKER = {}

# Translation table for reverse complement
fs = 'ACGTMRYKVHDBU'
ts = 'TGCAKYRMBDHVT'
ttable = str.maketrans(fs, ts)

def open_maybe_gzip(p, mode="rt"):
    return gzip.open(p, mode) if str(p).endswith(".gz") else open(p, mode)

# ----------------- Utilities -----------------
def reverse_complement(sequence):
    sequence = sequence.upper()
    return sequence[-1::-1].translate(ttable)

def dna2aa(sequence):
    gencode = {
        'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
        'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
        'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
        'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
        'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
        'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
        'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
        'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_','TGC':'C','TGT':'C','TGA':'_','TGG':'W',
    }
    return ''.join(
        [gencode[sequence[n:n+3]] for n in range(0, len(sequence), 3) if sequence[n:n+3] in gencode]
    )

def mismatches(peptide):
    mismatch = '[A-Z]'
    newpeps = []
    r = len(peptide)
    peptide2 = [letter for letter in peptide]
    for n in range(r):
        for m in range(r):
            peptide3 = peptide2[:]
            peptide3[n] = mismatch
            peptide3[m] = mismatch
            newpeps.append(''.join(peptide3))
    return '|'.join(set(newpeps))

def expand(peptide):
    if peptide.count('X') > 5:
        return []
    flag = 1
    oldlist = [peptide]
    length = len(peptide)
    aa = [l for l in 'ILVFMCAGPTSYWQNHEDKR']
    while flag == 1:
        newlist = []
        flag = 0
        for one in oldlist:
            for n in range(length):
                if one[n] == 'X':
                    flag = 1
                    x = one[n]
                    for l in aa:
                        newlist.append(one.replace(x, l, 1))
                    continue
        if newlist:
            oldlist = set(newlist)
    return oldlist

def levenshtein(arg):
    a, b, c = arg
    n, m = len(a), len(b)
    if n > m:
        a, b = b, a
        n, m = m, n
    current = list(range(n+1))
    conseq = 0
    for i in range(1, m+1):
        previous, current = current, [i] + [0]*n
        for j in range(1, n+1):
            if a[j-1] == b[i-1]:
                current[j] = previous[j-1]
            else:
                current[j] = min(previous[j] + 1, current[j-1] + 1, previous[j-1] + 1.9)
        for o in range(n):
            if previous[o] + 1 != current[o]:
                conseq = -1
                break
        if conseq == -1:
            conseq = 0
        else:
            conseq += 1
        f = min(current) - conseq
        if f > c:
            return f
    if conseq == 1:
        conseq = 0
    elif conseq > 1:
        conseq -= 2
    return current[n] - conseq

def levendata(one, dnafiles, etfile):
    """Merge related DNA variants (parent process)."""
    newdna = dnafiles
    etalons = {}
    mindna = 200
    dnadic = {}

    for bigpep in etfile:
        if bigpep in list(newdna.keys()):
            item = newdna[bigpep]
            for dna in item:
                if item[dna] > mindna:
                    etalons[dna] = item[dna]
                dnadic[dna] = [dna]

    peplist = dict((dna, pep) for pep in newdna for dna in newdna[pep])
    dnalist = dict((dna, newdna[pep][dna]) for pep in newdna for dna in newdna[pep])

    for item in sorted(list(dnalist.items()), key=itemgetter(1), reverse=True):
        dna = item[0]
        value = item[1]
        if dna in etalons:
            continue
        minlsd = args2.etaloncutoff + 1.0
        etargs = [(et, dna, minlsd) for et in etalons if etalons[et] > value * mindna]
        if args2.mp > 1 and pool is not None:
            minlsd_vals = pool.map(levenshtein, etargs)
            minlsdlist = [(minlsd_vals[n], etargs[n][0]) for n in range(len(etargs))]
        else:
            minlsdlist = [(levenshtein(et), et[0]) for et in etargs]
        for lsd, et in minlsdlist:
            if lsd < minlsd:
                minlsd = lsd
                minmatch = et
        if minlsd <= args2.etaloncutoff:
            dnadic[minmatch].append(dna)
            dnalist[dna] = float(dnalist[dna])
        else:
            dnadic[dna] = [dna]
            if item[1] > mindna and not args2.oldetalon:
                etalons[dna] = value

    newdna2 = {}
    if args2.etalonremove:
        for pep in dnafiles:
            for dna in dnafiles[pep]:
                if isinstance(dnalist[dna], float):
                    continue
                if pep not in newdna2:
                    newdna2[pep] = {}
                newdna2[pep][dna] = dnafiles[pep][dna]
    else:
        for dna in dnadic:
            pep = peplist[dna]
            if pep not in newdna2:
                newdna2[pep] = {}
            for dna2 in dnadic[dna]:
                newdna2[pep][dna2] = dnalist[dna2]

    two = defaultdict(int)
    for pep in newdna2:
        f = 0
        for dna in newdna2[pep]:
            f += newdna2[pep][dna]
        two[pep] = f
    return two, newdna2

# ----------------- Worker init -----------------
def _init_worker(cfg: dict):
    """Runs once per child process; builds args + regex inside worker."""
    a = types.SimpleNamespace(**cfg)

    def mismatches_local(peptide):
        mismatch = '[A-Z]'
        newpeps = []
        r = len(peptide)
        peptide2 = [letter for letter in peptide]
        for n in range(r):
            for m in range(r):
                peptide3 = peptide2[:]
                peptide3[n] = mismatch
                peptide3[m] = mismatch
                newpeps.append(''.join(peptide3))
        return '|'.join(set(newpeps))

    a.aa_form = re.compile(r'^([A-Z_]{0,'+str(a.stops)+'}[A-Z]*)')

    bb_form = False
    if a.his:
        bb_form = re.compile(r'^([A-Z_]{7})('+mismatches_local('ENLYFQR')+')[H]{5}[A-Z]*')
    elif a.cx4c:
        bb_form = re.compile(r'^(C[A-Z_]{4}C[A-Z]*)')
    elif a.cx7c:
        bb_form = re.compile(r'^(C[A-Z_]{7}C[A-Z]*)')
    elif a.cx8c:
        bb_form = re.compile(r'^(C[A-Z_]{8}C[A-Z]*)')
    elif a.cx10c:
        bb_form = re.compile(r'^(C[A-Z_]{10}C[A-Z]*)')
    elif not a.nocys:
        bb_form = re.compile(r'^(C[A-Z_]{3,}C[A-Z]*)')
    a.bb_form = bb_form

    a.mlg_form1 = re.compile(r'(MLGDPNS)([A-Z_]+)')
    a.mlg_form2 = re.compile(r'('+mismatches_local('MLGDPN')+')S([A-Z_]+)')

    _WORKER["args"] = a

# ----------------- Fallback for hot-reload/single-process -----------------
def _ensure_worker_args_from_globals():
    """Rebuild worker args from global args2 if worker state is empty."""
    if "args" in _WORKER and _WORKER["args"] is not None:
        return _WORKER["args"]

    a = types.SimpleNamespace(
        fasta=getattr(args2, "fasta", False),
        stops=getattr(args2, "stops", 0),
        nocys=getattr(args2, "nocys", False),
        rev=getattr(args2, "rev", False),
        fwd=getattr(args2, "fwd", False),
        limit=getattr(args2, "limit", 0),
        dnaspace=getattr(args2, "dnaspace", False),
        etaloncutoff=float(getattr(args2, "etaloncutoff", 4.0)),
        oldetalon=getattr(args2, "oldetalon", False),
        etalonremove=getattr(args2, "etalonremove", False),
        nnk=getattr(args2, "nnk", False),
        his=getattr(args2, "his", False),
        qual=getattr(args2, "qual", 0),
        cx7c=getattr(args2, "cx7c", False),
        cx8c=getattr(args2, "cx8c", False),
        cx10c=getattr(args2, "cx10c", False),
        cx4c=getattr(args2, "cx4c", False),
        dna=getattr(args2, "dna", False),
        etalon=bool(getattr(args2, "etalon", None)),
    )

    def mismatches_local(peptide):
        mismatch = '[A-Z]'
        newpeps = []
        r = len(peptide)
        peptide2 = [letter for letter in peptide]
        for n in range(r):
            for m in range(r):
                peptide3 = peptide2[:]
                peptide3[n] = mismatch
                peptide3[m] = mismatch
                newpeps.append(''.join(peptide3))
        return '|'.join(set(newpeps))

    a.aa_form = re.compile(r'^([A-Z_]{0,'+str(a.stops)+'}[A-Z]*)')

    bb_form = False
    if a.his:
        bb_form = re.compile(r'^([A-Z_]{7})('+mismatches_local('ENLYFQR')+')[H]{5}[A-Z]*')
    elif a.cx4c:
        bb_form = re.compile(r'^(C[A-Z_]{4}C[A-Z]*)')
    elif a.cx7c:
        bb_form = re.compile(r'^(C[A-Z_]{7}C[A-Z]*)')
    elif a.cx8c:
        bb_form = re.compile(r'^(C[A-Z_]{8}C[A-Z]*)')
    elif a.cx10c:
        bb_form = re.compile(r'^(C[A-Z_]{10}C[A-Z]*)')
    elif not a.nocys:
        bb_form = re.compile(r'^(C[A-Z_]{3,}C[A-Z]*)')
    a.bb_form = bb_form

    a.mlg_form1 = re.compile(r'(MLGDPNS)([A-Z_]+)')
    a.mlg_form2 = re.compile(r'('+mismatches_local('MLGDPN')+')S([A-Z_]+)')

    _WORKER["args"] = a
    return a

# ----------------- FASTQ/FASTA reader -----------------
def readfastagen(ffile):
    """Yields (seq, qual). FASTQ 4-line records, or FASTA concatenated lines."""
    name, seq = None, []
    n = 0
    qual = ''
    if getattr(args2, "fasta", False):
        for line in ffile:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield ''.join(seq), qual
                name, seq = line[1:], []
            else:
                seq.append(line)
    else:
        for line in ffile:
            line = line.rstrip()
            if n == 0:
                if name:
                    yield ''.join(seq), qual
                name, seq = line, []
            elif n == 1:
                seq.append(line)
            elif n == 3:
                qual = line
                n = -1
            n += 1
    if name:
        yield ''.join(seq), qual

# ----------------- Worker function -----------------
def process_line(arg):
    # arg: (seq, qual)
    line, q = arg
    a2 = _WORKER.get("args")
    if a2 is None:
        a2 = _ensure_worker_args_from_globals()

    pep = ''
    pepdna = ''
    cat = 0

    line = line.rstrip()
    line2 = reverse_complement(line)

    direction = []
    if a2.rev and a2.fwd:
        a2.rev = False
        a2.fwd = False
    if not a2.fwd:
        for p in [2, 1, 0]:
            direction.append([p, -1])
    if not a2.rev:
        for p in [2, 1, 0]:
            direction.append([p, 0])
    newdir = []
    for mlg in [1, 2]:
        for n in range(len(direction)):
            newdir.append(direction[n] + [mlg])
    direction = newdir

    for stoffset, order, mlg in direction:
        if order == -1:
            line3 = line2[stoffset:]
            qual = q[::-1][stoffset:]
        else:
            line3 = line[stoffset:]
            qual = q[stoffset:]

        peptide = dna2aa(line3)
        mlg_form = a2.mlg_form1 if mlg == 1 else a2.mlg_form2
        pepmatch = re.search(mlg_form, peptide)
        if pepmatch is None:
            continue
        cat = 1
        offset2 = (pepmatch.span()[0] * 3) + 21
        peptide2 = pepmatch.group(2)

        smallpep = re.search(a2.aa_form, peptide2)
        if smallpep is None:
            continue
        else:
            peptide3 = smallpep.group(1)
        cat = 2

        if a2.bb_form:
            smallpep = re.search(a2.bb_form, peptide3)
        if smallpep is None:
            continue
        else:
            spep = smallpep.group(1)
            cat = 3

        if a2.dna or getattr(a2, "etalon", False) or a2.nnk or a2.qual:
            dnastart = (smallpep.span()[0] * 3) + offset2
            dnaend2 = (smallpep.span()[1] * 3) + offset2
            dnaend = len(line3)
            nnkdnalen = dnaend2 - dnastart
            pepdna = line3[dnastart:dnaend]

            if a2.qual != 0:
                scores = [ord(n) - 33 for n in qual[dnastart:dnaend]]
                sc = [item for item in scores if item < a2.qual]
                if len(sc) > 3:
                    continue

            cat = 4
            if a2.nnk:
                error = 0
                for start in range(0, len(pepdna[:nnkdnalen]), 3):
                    triplet = pepdna[start:start+3]
                    if len(triplet) < 3:
                        break
                    if triplet[-1] in ('A', 'C'):
                        if triplet != "TGC":  # exception for Cys
                            error = 1
                if error == 1:
                    continue

        cat = 6 if mlg == 1 else 5
        pep = spep
        break

    if cat == 0 and len(line) < 70:
        cat = -1
    return pep, pepdna, cat

# ----------------- Output table -----------------
def sortoutput(files, dnafiles):
    peptides = {}
    keys = []
    lengths = []
    totalpeps = 0
    filerange = list(range(len(files)))
    for n in filerange:
        lengths.append(len(list(files[n].keys())))
        keys += list(files[n].keys())
        totalpeps += sum(files[n].values())
    keys = set(keys)
    uniqpeps = len(keys)
    maxlist = []
    if getattr(args2, "sort", False):
        for key in keys:
            peptides[key] = [files[n][key] for n in filerange]
            maxlist.append([key, [max(peptides[key])] + peptides[key]])
    else:
        for key in keys:
            peptides[key] = [files[n][key] for n in filerange]
            maxlist.append([key, peptides[key]])
    maxlist = sorted(maxlist, key=itemgetter(1), reverse=True)

    buffer = StringIO()
    print('"Unique peptides"'.ljust(18), end=' ', file=buffer)
    for item in lengths:
        print(repr(item).ljust(7), end=' ', file=buffer)
    print('\n', end=' ', file=buffer)
    print('"Total peptides"'.ljust(18), end=' ', file=buffer)
    for n in filerange:
        s = sum(files[n].values())
        print(repr(s).ljust(7), end=' ', file=buffer)
    print('\n', end=' ', file=buffer)
    print('\n"Total unique peptides:"', uniqpeps, file=buffer)
    print('"Sum total peptides:"', totalpeps, file=buffer)
    print('"Peptides"\t', end=' ', file=buffer)
    for new in titles:
        print('"'+ new + '"\t', end=' ', file=buffer)
    print('\n', end=' ', file=buffer)

    limit = len(maxlist) if args2.limit == 0 or args2.limit > len(maxlist) else args2.limit
    l = 0
    while l < limit:
        item = maxlist[l]
        l += 1
        key = item[0]
        print(key.ljust(18), end=' ', file=buffer)
        for n in filerange:
            print(repr(peptides[key][n]).ljust(7), end=' ', file=buffer)
        print('\n', end=' ', file=buffer)

        if getattr(args2, "dna", False):
            maxdna = ''
            for k in filerange:
                try:
                    item2 = dnafiles[k][key]
                except KeyError:
                    continue
                for dnakey in sorted(list(item2.items()), key=itemgetter(1), reverse=True):
                    seq = dnakey[0]
                    if maxdna == '':
                        maxdna = seq
                    g = 3 if getattr(args2, "dnaspace", False) else len(seq)
                    if not isinstance(dnakey[1], int):
                        sc = dnakey[1]
                        seq2 = dna2aa(seq)
                        seq3 = seq2
                        sc2 = levenshtein((maxdna, seq, args2.etaloncutoff + 1))
                    else:
                        seq3 = ''
                        sc = dnakey[1]
                        sc2 = ''
                    for m in range(0, len(seq), g):
                        print(seq[m:m+g], end=' ', file=buffer)
                    print(sc, seq3, sc2, file=buffer)
    return buffer.getvalue()

# ----------------- Public entry point -----------------
def run_on_files(
    filenames,
    *,
    fasta=False,
    stops=0,
    nocys=False,
    rev=False,
    fwd=False,
    limit=0,
    peplist_path=None,
    excludepl=False,
    dna=False,
    dnaspace=False,
    etalon_path=None,
    etaloncutoff=4.0,
    oldetalon=False,
    etalonremove=False,
    nnk=False,
    his=False,
    qual=0,
    cx7c=False,
    cx8c=False,
    cx10c=False,
    cx4c=False,
    report=False,          # kept for API symmetry
    mp_workers=1,
):
    """
    Returns: dict[str filename_base -> str output_text]
    """
    class A:  # argparse-like namespace
        pass
    global args2, titles, pool
    args2 = A()
    args2.sort = False
    args2.fasta = fasta
    args2.stops = stops
    args2.nocys = nocys
    args2.rev = rev
    args2.fwd = fwd
    args2.limit = limit
    args2.output = None
    args2.peplist = peplist_path is not None
    args2.excludepl = excludepl
    args2.dna = dna or dnaspace
    args2.dnaspace = dnaspace
    args2.etalon = (etalon_path if etalon_path else None)
    args2.etaloncutoff = float(etaloncutoff)
    args2.oldetalon = oldetalon
    args2.etalonremove = etalonremove
    args2.nnk = nnk
    args2.his = his
    args2.qual = qual
    args2.cx7c = cx7c
    args2.cx8c = cx8c
    args2.cx10c = cx10c
    args2.cx4c = cx4c
    args2.report = report
    args2.mp = mp_workers

    # Peplist
    plist = None
    if peplist_path:
        plist = []
        with open(peplist_path, 'r') as plfile:
            for line in plfile:
                if 'X' in line:
                    plist += expand(line.rstrip())
                else:
                    plist.append(line.rstrip())

    # Etalon list
    etlist = []
    if args2.etalon and not args2.oldetalon:
        with open(etalon_path, 'r') as etfile:
            for line in etfile:
                etlist.append(line.rstrip().split()[0])

    # Picklable worker config (booleans/numbers only)
    cfg = dict(
        fasta=args2.fasta,
        stops=args2.stops,
        nocys=args2.nocys,
        rev=args2.rev,
        fwd=args2.fwd,
        limit=args2.limit,
        dnaspace=args2.dnaspace,
        etaloncutoff=float(args2.etaloncutoff),
        oldetalon=args2.oldetalon,
        etalonremove=args2.etalonremove,
        nnk=args2.nnk,
        his=args2.his,
        qual=args2.qual,
        cx7c=args2.cx7c,
        cx8c=args2.cx8c,
        cx10c=args2.cx10c,
        cx4c=args2.cx4c,
        dna=args2.dna,
        etalon=bool(args2.etalon),
    )

    # Use real processes (no thread fallback)
    use_threads = False
    if args2.mp > 1:
        if use_threads:
            from multiprocessing.dummy import Pool as ThreadPool
            pool = ThreadPool(args2.mp)
            _init_worker(cfg)
        else:
            pool = mp.Pool(processes=args2.mp, initializer=_init_worker, initargs=(cfg,))
    else:
        pool = None
        _init_worker(cfg)  # seed for single-process path

    results_per_file = {}

    for filename in filenames:
        tmp = defaultdict(int)
        tmpdna = defaultdict(dict)
        cats = defaultdict(int)

        with open_maybe_gzip(filename, "rt") as fastq:
            iterable = readfastagen(fastq)
            if args2.mp > 1 and pool is not None and not use_threads:
                results_iter = pool.imap_unordered(process_line, iterable, chunksize=512)
            elif args2.mp > 1 and use_threads:
                results_iter = pool.imap(process_line, iterable, chunksize=512)
            else:
                results_iter = (process_line(item) for item in iterable)

            for spep, pepdna, cat in results_iter:
                if pepdna != '':
                    if spep in tmpdna:
                        tmpdna[spep][pepdna] += 1
                    else:
                        tmpdna[spep] = defaultdict(int)
                        tmpdna[spep][pepdna] += 1
                if spep != '':
                    tmp[spep] += 1
                cats[cat] += 1

        if args2.etalon and not args2.oldetalon:
            tmp, tmpdna = levendata(tmp, tmpdna, etlist)

        if peplist_path:
            tmp2 = tmp
            tmp = defaultdict(int)
            if args2.excludepl:
                for key in list(tmp2.keys()):
                    if key in plist:
                        pass
                    else:
                        tmp[key] = tmp2[key]
            else:
                for key in list(tmp2.keys()):
                    if key in plist:
                        tmp[key] = tmp2[key]

        files = [tmp]
        dnafiles = [tmpdna] if args2.dna else [{}]

        # Per-file title (override global for this output)
        base = path.basename(filename).rsplit('.', 1)[0]
        globals()['titles'] = [base]

        single_out = sortoutput(files, dnafiles)
        results_per_file[path.basename(filename).rsplit('.', 1)[0]] = single_out

    if pool:
        pool.close()

    return results_per_file
