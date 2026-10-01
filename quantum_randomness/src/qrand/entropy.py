"""T3 — Entropia a blocchi H_n, entropia condizionale h_n e G-test di ordine di Markov."""

import numpy as np

from .bits import iter_block_codes, permuted
from .stats import g_test_table, miller_madow_entropy

CAP = 2 ** 27
N_MAX = 20
MARKOV_K = (1, 2, 3, 4, 6, 8, 10, 12, 14, 16)
N_SURROGATES = 5


def block_entropies(bits, n_max=N_MAX, markov_k=()):
    """H_n (Miller–Madow) per n = 1..n_max+1 e tabelle (contesto x bit successivo) per i k richiesti."""
    H, tables = {}, {}
    n_top = max(n_max + 1, max(markov_k, default=0) + 1)
    for n, codes in iter_block_codes(bits, n_top):
        counts = np.bincount(codes, minlength=2 ** n)
        H[n] = miller_madow_entropy(counts)
        if n - 1 in markov_k:
            tables[n - 1] = counts.reshape(2 ** (n - 1), 2)
    return H, tables


def conditional(H, n_max=N_MAX):
    """h_0 = H_1, h_n = H_{n+1} - H_n per n = 1..n_max."""
    h = [H[1]] + [H[n + 1] - H[n] for n in range(1, n_max + 1)]
    return np.array(h)


def run(bits, label=""):
    b = bits[:CAP]
    N = len(b)
    p = float(b.mean())
    q = min(p, 1 - p)
    ks = [k for k in MARKOV_K if q > 0 and N * q ** (k + 1) >= 5]
    n_max = min(N_MAX, max(1, len(b) - 2))
    H, tables = block_entropies(b, n_max, ks)
    pvals = []
    for k in ks:
        g, dof, pv = g_test_table(tables[k])
        pvals.append({"test": "T3", "name": f"markov_k{k}", "p": pv, "stat": g, "dof": dof})
    h = conditional(H, n_max - 1)
    band = [conditional(block_entropies(permuted(b), n_max)[0], n_max - 1)
            for _ in range(N_SURROGATES)]
    reliable = [2 ** (n + 1) <= N / 10 for n in range(len(h))]
    desc = {"N": N, "p_hat": p, "H": [H[n] for n in range(1, n_max + 1)], "h": h.tolist(),
            "h_surrogati": np.array(band).tolist(), "affidabile": reliable,
            "markov_k_testati": ks}
    return pvals, desc
