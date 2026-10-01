"""T2 — Stima di min-entropia NIST SP 800-90B (tool ufficiale `ea_non_iid`, modalità non-IID)."""

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from . import TOOLS_DIR
from .bits import bernoulli_bits

EA_NON_IID = TOOLS_DIR / "SP800-90B_EntropyAssessment" / "cpp" / "ea_non_iid"
N_SAMPLES = 1_000_000
N_REF = 19


def available():
    return EA_NON_IID.exists()


def min_entropy(bits):
    """(H_min complessiva, stime per estimatore) dei primi 10^6 bit, un simbolo per byte."""
    with tempfile.TemporaryDirectory(prefix="ea_") as tmp:
        sym = Path(tmp) / "data.bin"
        out = Path(tmp) / "out.json"
        np.asarray(bits[:N_SAMPLES], dtype=np.uint8).tofile(sym)
        subprocess.run([str(EA_NON_IID), "-q", "-i", "-o", str(out), str(sym), "1"],
                       cwd=tmp, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        report = json.loads(out.read_text())
    return _extract_hmin(report)


def _extract_hmin(report):
    """H_min complessiva (caso "Overall") e stime dei singoli estimatori."""
    per_test, overall = {}, None
    for case in report["testCases"]:
        if case["testCaseDesc"] == "Overall":
            overall = float(case["hAssessed"])
        else:
            val = case.get("hOriginal", case.get("tTupleRes", case.get("lrsRes")))
            per_test[case["testCaseDesc"]] = float(val)
    if overall is None:
        raise KeyError("caso 'Overall' assente nel JSON di ea_non_iid")
    return overall, per_test


def run(bits, label=""):
    """Descrittivo: H_min dei dati e rango rispetto a 19 surrogati i.i.d. Bernoulli(p̂)."""
    if len(bits) < N_SAMPLES:
        return [], {"applicabile": False, "motivo": f"N = {len(bits)} < 10^6"}
    p_hat = float(np.mean(bits[:N_SAMPLES]))
    h, per_test = min_entropy(bits)
    refs = [min_entropy(bernoulli_bits(N_SAMPLES, p_hat))[0] for _ in range(N_REF)]
    p_rank = (1 + sum(r <= h for r in refs)) / (N_REF + 1)
    return [], {
        "applicabile": True, "p_hat": p_hat, "H_min": h, "stimatori": per_test, "H_ref": refs,
        "H_ref_min": min(refs), "H_teorica": float(-np.log2(max(p_hat, 1 - p_hat))),
        "p_rank": p_rank, "segnalata": bool(h < min(refs)),
    }
