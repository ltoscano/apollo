"""Strumenti statistici comuni: correzioni per test multipli, combinazioni, G-test."""

import numpy as np
from scipy import stats


def benjamini_hochberg(pvals, q):
    """Restituisce (p aggiustati BH, maschera dei rifiuti a livello q)."""
    p = np.asarray(pvals, dtype=float)
    m = p.size
    if m == 0:
        return p, np.zeros(0, dtype=bool)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.minimum(adj_sorted, 1.0)
    return adj, adj <= q


def bonferroni(pvals, alpha):
    p = np.asarray(pvals, dtype=float)
    adj = np.minimum(p * p.size, 1.0)
    return adj, adj <= alpha


def fisher_combine(pvals):
    p = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1.0)
    if p.size == 0:
        return np.nan
    return float(stats.chi2.sf(-2 * np.log(p).sum(), 2 * p.size))


def g_test_table(table):
    """G-test di indipendenza su una tabella di contingenza 2-D (righe/colonne vuote rimosse)."""
    t = np.asarray(table, dtype=float)
    t = t[t.sum(axis=1) > 0][:, t.sum(axis=0) > 0]
    if t.shape[0] < 2 or t.shape[1] < 2:
        return 0.0, 0, 1.0
    n = t.sum()
    expected = np.outer(t.sum(axis=1), t.sum(axis=0)) / n
    nz = t > 0
    g = 2.0 * float((t[nz] * np.log(t[nz] / expected[nz])).sum())
    dof = (t.shape[0] - 1) * (t.shape[1] - 1)
    return g, dof, float(stats.chi2.sf(g, dof))


def chi2_gof(observed, expected_prob, min_expected=5.0, n_fitted=0):
    """Chi-quadro di bontà di adattamento; le celle con atteso basso vengono fuse in coda
    (dopo ordinamento per atteso decrescente)."""
    obs = np.asarray(observed, dtype=float)
    prob = np.asarray(expected_prob, dtype=float)
    n = obs.sum()
    exp = prob / prob.sum() * n
    order = np.argsort(-exp)
    obs, exp = obs[order], exp[order]
    keep = exp >= min_expected
    if keep.all():
        o, e = obs, exp
    else:
        o = np.append(obs[keep], obs[~keep].sum())
        e = np.append(exp[keep], exp[~keep].sum())
        if e[-1] < min_expected and len(e) > 1:
            o[-2] += o[-1]
            e[-2] += e[-1]
            o, e = o[:-1], e[:-1]
    dof = len(o) - 1 - n_fitted
    if dof < 1:
        return np.nan, 0, np.nan
    chi2 = float(((o - e) ** 2 / e).sum())
    return chi2, dof, float(stats.chi2.sf(chi2, dof))


def surrogate_t_pvalue(value, surrogates, alternative="less"):
    """p-value di predizione t di `value` rispetto a M surrogati (t con M-1 gradi di libertà)."""
    s = np.asarray(surrogates, dtype=float)
    m = s.size
    sd = s.std(ddof=1)
    if sd == 0:
        sd = 1e-12
    t = (value - s.mean()) / (sd * np.sqrt(1 + 1 / m))
    if alternative == "less":
        return float(t), float(stats.t.cdf(t, m - 1))
    if alternative == "greater":
        return float(t), float(stats.t.sf(t, m - 1))
    return float(t), float(2 * stats.t.sf(abs(t), m - 1))


def miller_madow_entropy(counts):
    """Entropia (bit) plug-in con correzione di Miller–Madow."""
    c = np.asarray(counts, dtype=float)
    c = c[c > 0]
    n = c.sum()
    if n == 0:
        return 0.0
    p = c / n
    return float(-(p * np.log2(p)).sum() + (c.size - 1) / (2 * n * np.log(2)))
