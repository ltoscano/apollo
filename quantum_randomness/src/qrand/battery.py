"""Esecuzione della batteria pre-registrata su una sequenza e correzione per test multipli."""

import time

import numpy as np

from . import ALPHA, autocorr, compress, entropy, linguistic, nist22, sp90b
from .stats import benjamini_hochberg, bonferroni

BIT_TESTS = {
    "T1": nist22.run,
    "T2": sp90b.run,
    "T3": entropy.run,
    "T4": autocorr.run,
    "T5": compress.run,
    "T6": linguistic.run,
}


def run_bits(bits, label, fair_null, tests=tuple(BIT_TESTS), log=print):
    """Esegue i test su una sequenza di bit. Con fair_null=False i p-value di T1 (SP 800-22,
    che assume Bernoulli(1/2)) sono marcati come descrittivi e restano fuori dalla famiglia."""
    bits = np.ascontiguousarray(bits, dtype=np.uint8)
    out = {"label": label, "N": int(len(bits)), "p_hat": float(bits.mean()) if len(bits) else None,
           "fair_null": fair_null, "pvalues": [], "descrittivi": {}}
    for name in tests:
        if name == "T1" and not nist22.available():
            out["descrittivi"]["T1"] = {"saltato": "sts-2.1.2 non compilato"}
            continue
        if name == "T2" and not sp90b.available():
            out["descrittivi"]["T2"] = {"saltato": "ea_non_iid non compilato"}
            continue
        t0 = time.time()
        pv, desc = BIT_TESTS[name](bits, label)
        confirmatory = not (name == "T1" and not fair_null)
        for rec in pv:
            rec["sequenza"] = label
            rec["confermativo"] = confirmatory
        out["pvalues"] += pv
        out["descrittivi"][name] = desc
        log(f"  [{label}] {name}: {len(pv)} p-value, {time.time() - t0:.1f}s")
    return out


def run_symbols(symbols, alphabet, label, log=print):
    """Sequenze di simboli (`joint`, Sycamore `hw`), nullo "prove i.i.d.": G-test di Markov sui
    simboli, T4 sui valori dei simboli centrati, T5 con un byte per simbolo e surrogati per
    permutazione dei simboli (deviazione 1 della pre-registrazione)."""
    s = np.asarray(symbols, dtype=np.int64)
    out = {"label": label, "N": int(s.size), "alfabeto": alphabet, "fair_null": False,
           "pvalues": [], "descrittivi": {}}
    pv = linguistic.markov_symbols(s, alphabet)
    p, d = autocorr.run(s, label)
    pv += p
    out["descrittivi"]["T4"] = d
    p, d = compress.run(s, label, symbols=True)
    pv += p
    out["descrittivi"]["T5"] = d
    for rec in pv:
        rec["sequenza"] = label
        rec["confermativo"] = True
    out["pvalues"] = pv
    log(f"  [{label}] simboli: {len(pv)} p-value")
    return out


def correct(results, q=ALPHA, alpha=ALPHA):
    """Applica BH e Bonferroni alla famiglia dei p-value confermativi di un dataset."""
    records = [r for res in results for r in res["pvalues"] if r.get("confermativo")]
    p = np.array([r["p"] for r in records], dtype=float)
    p = np.where(np.isnan(p), 1.0, p)
    bh_adj, bh_rej = benjamini_hochberg(p, q)
    bf_adj, bf_rej = bonferroni(p, alpha)
    for r, a, ra, b, rb in zip(records, bh_adj, bh_rej, bf_adj, bf_rej):
        r.update(p_bh=float(a), sig_bh=bool(ra), p_bonf=float(b), sig_bonf=bool(rb))
    return {"m": len(records), "scoperte_bh": int(bh_rej.sum()), "scoperte_bonf": int(bf_rej.sum()),
            "records": records}
