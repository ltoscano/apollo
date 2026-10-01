"""T6 — Analisi "linguistica" sul modello di Doyle & McCowan: parole a lunghezza fissa
e variabile, legge di Zipf, entropie condizionali di ordine superiore."""

import numpy as np

from .bits import permuted, words
from .stats import chi2_gof, g_test_table, miller_madow_entropy, surrogate_t_pvalue

CAP = 2 ** 27
N_SURROGATES = 10
ZIPF_RANKS = 1000
COND_ORDERS = 3
MAX_TOKEN_BINS = 16


def rank_frequency(tokens):
    counts = np.bincount(tokens) if tokens.size else np.zeros(0, int)
    counts = np.sort(counts[counts > 0])[::-1]
    return counts


def zipf_slope(counts, max_rank=ZIPF_RANKS):
    c = counts[:max_rank]
    if c.size < 3:
        return np.nan
    r = np.arange(1, c.size + 1)
    return float(np.polyfit(np.log(r), np.log(c), 1)[0])


def word_probs(w, p):
    ones = np.array([bin(i).count("1") for i in range(2 ** w)])
    return p ** ones * (1 - p) ** (w - ones)


def conditional_entropies(tokens, orders=COND_ORDERS):
    """h_k = H(k+1-grammi) - H(k-grammi) dei token, con Miller–Madow; k = 0..orders."""
    H = [0.0]
    for n in range(1, orders + 2):
        if tokens.size < n:
            break
        grams = np.stack([tokens[i:tokens.size - n + 1 + i] for i in range(n)], axis=1)
        _, counts = np.unique(grams, axis=0, return_counts=True)
        H.append(miller_madow_entropy(counts))
    return [H[i + 1] - H[i] for i in range(len(H) - 1)]


def token_lengths(bits, sep):
    pos = np.flatnonzero(bits == sep)
    return np.diff(pos)


def _length_bins(n_tokens, q, joint):
    """Numero di classi B (1..B-1 e coda >= B) tale che gli attesi siano >= 5
    (sulle marginali o, se joint, sulle celle della tabella congiunta)."""
    best = 0
    for B in range(2, MAX_TOKEN_BINS + 1):
        probs = np.append(q * (1 - q) ** np.arange(B - 1), (1 - q) ** (B - 1))
        need = n_tokens * probs.min() ** 2 if joint else n_tokens * probs.min()
        if need >= 5:
            best = B
    return best


def fixed_words(bits, w, surrogate_bits):
    tok = words(bits, w)
    p = float(bits.mean())
    pvals, desc = [], {"n_parole": int(tok.size)}
    if tok.size >= 5 * 2 ** w:
        counts = np.bincount(tok, minlength=2 ** w)
        chi2, dof, pv = chi2_gof(counts, word_probs(w, p), n_fitted=1)
        pvals.append({"test": "T6", "name": f"parole{w}_freq", "p": pv, "stat": chi2, "dof": dof})
    rf = rank_frequency(tok)
    slope = zipf_slope(rf)
    s_rf = [rank_frequency(words(s, w)) for s in surrogate_bits]
    s_slopes = [zipf_slope(c) for c in s_rf]
    if np.isfinite(slope) and np.all(np.isfinite(s_slopes)):
        t, pv = surrogate_t_pvalue(slope, s_slopes, "two-sided")
        pvals.append({"test": "T6", "name": f"parole{w}_zipf", "p": pv, "stat": t})
    desc.update(zipf_pendenza=slope, zipf_pendenze_surrogati=s_slopes,
                rank_freq=rf[:ZIPF_RANKS].tolist(), rank_freq_surrogato=s_rf[0][:ZIPF_RANKS].tolist())
    if w == 8:
        if tok.size >= 20 * 2 ** 16:
            table = np.bincount(tok[:-1] * 256 + tok[1:], minlength=2 ** 16).reshape(256, 256)
            g, dof, pv = g_test_table(table)
            pvals.append({"test": "T6", "name": "parole8_ordine1", "p": pv, "stat": g, "dof": dof})
        sub = tok[: 2 ** 22]
        desc["h_cond"] = conditional_entropies(sub)
        desc["h_cond_surrogati"] = [conditional_entropies(words(s, 8)[: 2 ** 22])
                                    for s in surrogate_bits[:3]]
    return pvals, desc


def variable_tokens(bits, sep):
    lengths = token_lengths(bits, sep)
    q = float((bits == sep).mean())
    pvals, desc = [], {"n_token": int(lengths.size), "q": q}
    if lengths.size < 10 or q in (0.0, 1.0):
        return pvals, desc
    B = _length_bins(lengths.size, q, joint=False)
    if B >= 3:
        obs = np.bincount(np.minimum(lengths, B), minlength=B + 1)[1:]
        probs = np.append(q * (1 - q) ** np.arange(B - 1), (1 - q) ** (B - 1))
        chi2, dof, pv = chi2_gof(obs, probs, n_fitted=1)
        pvals.append({"test": "T6", "name": f"token_sep{sep}_geom", "p": pv, "stat": chi2,
                      "dof": dof})
    Bj = _length_bins(lengths.size, q, joint=True)
    if Bj >= 2:
        capped = np.minimum(lengths, Bj) - 1
        table = np.bincount(capped[:-1] * Bj + capped[1:], minlength=Bj * Bj).reshape(Bj, Bj)
        g, dof, pv = g_test_table(table)
        pvals.append({"test": "T6", "name": f"token_sep{sep}_ordine1", "p": pv, "stat": g,
                      "dof": dof})
    desc["rank_freq"] = rank_frequency(lengths)[:ZIPF_RANKS].tolist()
    return pvals, desc


def run(bits, label=""):
    b = bits[:CAP]
    surrogates = [permuted(b) for _ in range(N_SURROGATES)]
    pvals, desc = [], {"N": len(b)}
    for w in (8, 16):
        pv, d = fixed_words(b, w, surrogates)
        pvals += pv
        desc[f"parole{w}"] = d
    for sep in (1, 0):
        pv, d = variable_tokens(b, sep)
        pvals += pv
        desc[f"token_sep{sep}"] = d
    return pvals, desc


def markov_symbols(symbols, alphabet, orders=(1, 2)):
    """G-test di ordine di Markov su una sequenza di simboli (per `joint` e Sycamore `hw`)."""
    s = np.asarray(symbols, dtype=np.int64)
    pvals = []
    for k in orders:
        n = s.size - k
        if n < 5 * alphabet ** (k + 1):
            continue
        ctx = np.zeros(n, dtype=np.int64)
        for i in range(k):
            ctx = ctx * alphabet + s[i:i + n]
        table = np.bincount(ctx * alphabet + s[k:], minlength=alphabet ** (k + 1))
        g, dof, pv = g_test_table(table.reshape(alphabet ** k, alphabet))
        pvals.append({"test": "T3", "name": f"markov_simboli_k{k}", "p": pv, "stat": g, "dof": dof})
    return pvals
