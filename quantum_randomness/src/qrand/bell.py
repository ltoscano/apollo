"""Sanity check fisico dei test di Bell: CHSH (Delft) e CH/Eberhard con test martingala (NIST)."""

import itertools

import numpy as np
from scipy import stats

DELFT_PUBLISHED = {"n": 245, "k": 196, "p": 0.039, "S": 2.42, "S_err": 0.20}
NIST_P_THRESHOLD = 1e-6  # stesso ordine del p <= 2.3e-7 dichiarato (pre-registrazione, §4)


def correlators(a_set, b_set, a_out, b_out):
    E, n = np.zeros((2, 2)), np.zeros((2, 2), dtype=int)
    for x, y in itertools.product((0, 1), repeat=2):
        sel = (a_set == x) & (b_set == y)
        n[x, y] = sel.sum()
        if n[x, y]:
            E[x, y] = np.mean(np.where(a_out[sel] == b_out[sel], 1.0, -1.0))
    return E, n


def chsh_all_conventions(a_set, b_set, a_out, b_out):
    """S, errore standard, vittorie k e p-value binomiale P(K >= k | Bin(n, 3/4)) per le 8
    convenzioni: correlatore con segno meno in (x0, y0) e segno globale +/-."""
    a_set, b_set, a_out, b_out = (np.asarray(v, dtype=int) for v in (a_set, b_set, a_out, b_out))
    E, n = correlators(a_set, b_set, a_out, b_out)
    se_terms = np.where(n > 0, (1 - E ** 2) / np.maximum(n, 1), 0.0)
    out = []
    for (x0, y0), sign in itertools.product(itertools.product((0, 1), repeat=2), (1, -1)):
        signs = np.ones((2, 2))
        signs[x0, y0] = -1
        S = sign * float((signs * E).sum())
        minus = (a_set == x0) & (b_set == y0)
        target = np.where(minus, 1, 0) if sign == 1 else np.where(minus, 0, 1)
        k = int(((a_out ^ b_out) == target).sum())
        N = len(a_out)
        out.append({"meno_su": [x0, y0], "segno": sign, "S": S,
                    "S_err": float(np.sqrt(se_terms.sum())), "n": N, "k": k,
                    "p": float(stats.binom.sf(k - 1, N, 0.75))})
    return out, E, n


def delft_sanity(a_set, b_set, a_out, b_out, ref=DELFT_PUBLISHED):
    conv, E, n = chsh_all_conventions(a_set, b_set, a_out, b_out)
    ok = [c for c in conv
          if c["n"] == ref["n"] and c["k"] == ref["k"] and abs(c["p"] - ref["p"]) < 5e-4
          and abs(c["S"] - ref["S"]) < 0.01]
    return {"passato": len(ok) == 1, "convenzioni": conv, "corrispondenti": ok,
            "E": E.tolist(), "n_per_setting": n.tolist(), "pubblicato": ref}


def ch_eberhard_J(a_set, b_set, a_out, b_out):
    """J = P(11|00) - P(10|01) - P(01|10) - P(11|11); J > 0 viola il realismo locale."""
    a_set, b_set, a_out, b_out = (np.asarray(v, dtype=int) for v in (a_set, b_set, a_out, b_out))

    def P(x, y, a, b):
        sel = (a_set == x) & (b_set == y)
        return float(np.mean((a_out[sel] == a) & (b_out[sel] == b))) if sel.any() else np.nan

    return P(0, 0, 1, 1) - P(0, 1, 1, 0) - P(1, 0, 0, 1) - P(1, 1, 1, 1)


def bell_function(a_set, b_set, a_out, b_out):
    """B per prova, con E_LR[B] <= 0 se i setting sono uniformi e indipendenti (p = 1/4)."""
    x, y, a, b = a_set, b_set, a_out, b_out
    return 4.0 * (((x == 0) & (y == 0) & (a == 1) & (b == 1)).astype(float)
                  - ((x == 0) & (y == 1) & (a == 1) & (b == 0))
                  - ((x == 1) & (y == 0) & (a == 0) & (b == 1))
                  - ((x == 1) & (y == 1) & (a == 1) & (b == 1)))


def martingale_pvalue(a_set, b_set, a_out, b_out, train_frac=0.1):
    """p-value di un test supermartingala prod(1 + lam*B): lam ottimizzato sul primo 10% delle
    prove (training) e applicato al restante 90%. Vale p <= 1/prod (disuguaglianza di Ville)
    sotto realismo locale con setting uniformi e indipendenti."""
    a_set, b_set, a_out, b_out = (np.asarray(v, dtype=int) for v in (a_set, b_set, a_out, b_out))
    B = bell_function(a_set, b_set, a_out, b_out)
    m = int(len(B) * train_frac)
    train, test = B[:m], B[m:]
    vals, cnt = np.unique(train, return_counts=True)
    grid = np.linspace(0, 0.25, 2001)[1:-1]  # 1 + lam*B >= 0 richiede lam <= 1/4
    logs = [(cnt * np.log1p(l * vals)).sum() for l in grid]
    lam = float(grid[int(np.argmax(logs))])
    vt, ct = np.unique(test, return_counts=True)
    log_t = float((ct * np.log1p(lam * vt)).sum())
    return {"lambda": lam, "log_T": log_t, "p": float(min(1.0, np.exp(-log_t))),
            "n_train": m, "n_test": int(len(test))}


def nist_sanity(a_set, b_set, a_out, b_out):
    J = ch_eberhard_J(a_set, b_set, a_out, b_out)
    mart = martingale_pvalue(a_set, b_set, a_out, b_out)
    return {"passato": bool(J > 0 and mart["p"] <= NIST_P_THRESHOLD), "J": J,
            "martingala": mart, "soglia_p": NIST_P_THRESHOLD,
            "pubblicato": {"p_min": 5.9e-9, "p_aggiustato": 2.3e-7}}


def simulate_qm_trials(n, rng, visibility=1.0, eff=1.0):
    """Prove sintetiche CHSH con statistica quantistica (stato di singoletto, angoli ottimali)
    — solo per i test di unità, mai come dato."""
    x = rng.integers(0, 2, n)
    y = rng.integers(0, 2, n)
    ang_a = np.where(x == 0, 0.0, np.pi / 2)
    ang_b = np.where(y == 0, np.pi / 4, -np.pi / 4)
    corr = visibility * np.cos(ang_a - ang_b)  # E = cos(a - b), segno meno su (1,1)
    same = rng.random(n) < (1 + corr) / 2
    a = rng.integers(0, 2, n)
    b = np.where(same, a, 1 - a)
    if eff < 1:
        a = a & (rng.random(n) < eff)
        b = b & (rng.random(n) < eff)
    return x, y, a, b
