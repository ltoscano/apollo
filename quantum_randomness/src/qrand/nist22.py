"""T1 — NIST SP 800-22 tramite la suite C ufficiale sts-2.1.2 (`assess`).

`assess` è interattivo e scrive i risultati in percorsi relativi alla directory corrente
(experiments/AlgorithmTesting/..., templates/...). Ogni esecuzione avviene quindi in una
directory temporanea con la struttura attesa; le risposte al menu vengono passate su stdin.
"""

import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import special

from . import TOOLS_DIR
from .bits import to_bytes
from .stats import fisher_combine

STS_DIR = TOOLS_DIR / "sts"
STREAM_LEN = 1_000_000
MAX_STREAMS = 100

TESTS = [
    "Frequency", "BlockFrequency", "CumulativeSums", "Runs", "LongestRun", "Rank", "FFT",
    "NonOverlappingTemplate", "OverlappingTemplate", "Universal", "ApproximateEntropy",
    "RandomExcursions", "RandomExcursionsVariant", "Serial", "LinearComplexity",
]
# Lunghezza minima raccomandata (SP 800-22 rev1a, §2); per NonOverlappingTemplate la soglia
# è una scelta conservativa (8 blocchi con abbastanza occorrenze attese del template).
MIN_LEN = {
    "Frequency": 100, "BlockFrequency": 100, "CumulativeSums": 100, "Runs": 100,
    "LongestRun": 128, "Rank": 38_912, "FFT": 1_000, "NonOverlappingTemplate": 100_000,
    "OverlappingTemplate": 1_000_000, "Universal": 387_840, "ApproximateEntropy": 100,
    "RandomExcursions": 1_000_000, "RandomExcursionsVariant": 1_000_000, "Serial": 100,
    "LinearComplexity": 1_000_000,
}
# Test che producono più p-value per sequenza: numero di file dataK.txt.
MULTI = {"CumulativeSums": 2, "Serial": 2, "NonOverlappingTemplate": 148,
         "RandomExcursions": 8, "RandomExcursionsVariant": 18}


def available():
    return (STS_DIR / "assess").exists()


def _prepare_workdir():
    work = Path(tempfile.mkdtemp(prefix="sts_"))
    (work / "templates").symlink_to(STS_DIR / "templates")
    for t in TESTS:
        (work / "experiments" / "AlgorithmTesting" / t).mkdir(parents=True)
    return work


def _read_floats(path):
    return np.array([float(x) for x in path.read_text().split()], dtype=float)


def _valid_excursion_mask(stats_path, header, k):
    """Sequenze per cui il test di escursioni casuali è applicabile (J >= 500). Se assess ha
    interrotto il test su qualche sequenza (troppi cicli: J > max(1000, n/100)), non scrive
    nulla per quella sequenza e la corrispondenza riga-sequenza è persa: si restituisce None."""
    blocks = stats_path.read_text().split(header)[1:]
    if len(blocks) != k:
        return None
    return np.array(["INSUFFICIENT" not in b for b in blocks])


def uniformity_pvalue(p):
    """P-value di uniformità NIST (10 classi, chi^2 con 9 g.d.l.), con atteso non arrotondato."""
    hist = np.bincount(np.minimum((np.asarray(p) * 10).astype(int), 9), minlength=10)
    exp = len(p) / 10
    chi2 = ((hist - exp) ** 2 / exp).sum()
    return float(special.gammaincc(4.5, chi2 / 2))


PARAM_TESTS = ["BlockFrequency", "NonOverlappingTemplate", "OverlappingTemplate",
               "ApproximateEntropy", "Serial", "LinearComplexity"]


def _assess_once(bits, n, k, selected, timeout):
    work = _prepare_workdir()
    try:
        data = work / "input.bin"
        # Il lettore binario di assess consuma parole di 4 byte: si aggiunge un riempimento,
        # che non viene letto (assess usa solo n*k bit).
        raw = to_bytes(bits[: n * k])
        data.write_bytes(raw + bytes(4 + (-len(raw)) % 4))
        if len(selected) == len(TESTS):
            choice = "1\n"
        else:
            choice = "0\n" + "".join("1" if t in selected else "0" for t in TESTS) + "\n"
        params = ""
        if n < STREAM_LEN:
            log2n = int(math.floor(math.log2(n)))
            menu = [t for t in PARAM_TESTS if t in selected]
            new = {"ApproximateEntropy": max(1, min(10, log2n - 6)),
                   "Serial": max(2, min(16, log2n - 3))}
            for t, value in new.items():
                if t in menu:
                    params += f"{menu.index(t) + 1}\n{value}\n"
        # Il menu dei parametri compare solo se è selezionato almeno un test parametrico.
        menu_end = "0\n" if any(t in PARAM_TESTS for t in selected) else ""
        answers = f"0\n{data}\n{choice}{params}{menu_end}{k}\n1\n"
        # assess termina con codice 1 anche quando va a buon fine e, se un p-value non è un numero
        # in [0, 1] (succede con sequenze corte), va in segfault nel rapporto finale, che qui non
        # si usa. Si verificano quindi i file per test.
        proc = subprocess.run([str(STS_DIR / "assess"), str(n)], input=answers.encode(), cwd=work,
                              capture_output=True, timeout=timeout)
        out, aborted = {}, []
        base = work / "experiments" / "AlgorithmTesting"
        for t in selected:
            d = base / t
            try:
                if t in MULTI:
                    cols = [_read_floats(d / f"data{i + 1}.txt") for i in range(MULTI[t])]
                    if t == "Serial" and k == 1 and not cols[0].size:
                        cols = [_read_floats(d / "results.txt")[i:i + 1] for i in range(2)]
                    mat = np.stack(cols, axis=1)
                else:
                    mat = _read_floats(d / "results.txt")[:, None]
            except (OSError, ValueError):
                mat = np.zeros((0, 1))
            if mat.shape[0] != k:
                if len(selected) > 1:
                    raise RuntimeError(f"assess (codice {proc.returncode}): {t} ha {mat.shape[0]} "
                                       f"righe invece di {k}")
                aborted.append(t)
                continue
            mat[(mat < 0) | (mat > 1)] = np.nan
            if t in ("RandomExcursions", "RandomExcursionsVariant"):
                header = "RANDOM EXCURSIONS VARIANT TEST" if t.endswith("Variant") \
                    else "RANDOM EXCURSIONS TEST"
                mask = _valid_excursion_mask(d / "stats.txt", header, k)
                if mask is None:
                    aborted.append(t)
                    mat[:] = np.nan
                else:
                    mat[~mask] = np.nan
            out[t] = mat
        return out, aborted
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_assess(bits, n, k, timeout=36000):
    """Esegue assess su k sequenze di n bit; restituisce ({test: matrice (k, n_pvalue)} con NaN
    dove il test non è applicabile, lista dei test interrotti). Con una sola sequenza corta ogni
    test gira in un processo separato, così un crash di assess non perde gli altri risultati."""
    selected = [t for t in TESTS if n >= MIN_LEN[t]]
    if not selected:
        return {}, []
    if k > 1:
        return _assess_once(bits, n, k, selected, timeout)
    out, aborted = {}, []
    for t in selected:
        o, a = _assess_once(bits, n, k, [t], timeout)
        out.update(o)
        aborted += a
    return out, aborted


def run(bits, label=""):
    """Esegue T1 secondo la pre-registrazione; restituisce (lista di p-value, descrittivi)."""
    N = len(bits)
    if N >= STREAM_LEN:
        n, k = STREAM_LEN, min(N // STREAM_LEN, MAX_STREAMS)
    else:
        n, k = N, 1
    res, aborted = run_assess(bits, n, k)
    pvals, desc = [], {"n": n, "k": k, "rows": [], "interrotti": aborted}
    for t, mat in res.items():
        if n < MIN_LEN[t]:
            continue
        for j in range(mat.shape[1]):
            col = mat[:, j]
            col = col[~np.isnan(col)]
            if col.size == 0:
                continue
            if k == 1:
                p, method = float(col[0]), "primo livello"
            elif col.size >= 55:
                p, method = uniformity_pvalue(col), "uniformità"
            else:
                p, method = fisher_combine(col), "Fisher"
            passed = int((col >= 0.01).sum())
            name = f"{t}#{j + 1}" if mat.shape[1] > 1 else t
            pvals.append({"test": "T1", "name": name, "p": p, "method": method})
            desc["rows"].append({"name": name, "p": p, "method": method,
                                 "pass": passed, "valid": int(col.size)})
    return pvals, desc
