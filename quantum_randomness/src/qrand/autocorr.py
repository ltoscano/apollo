"""T4 — Autocorrelazione (lag 1..10^4) e spettro di potenza (FFT)."""

import numpy as np
from scipy import stats

CAP = 2 ** 25
MAX_LAG = 10_000


def acf(x, max_lag):
    n = len(x)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, size)
    r = np.fft.irfft(f * np.conj(f), size)[: max_lag + 1]
    return r[1:] / r[0]


def fisher_g_pvalue(periodogram):
    """Test g di Fisher sul picco massimo; approssimazione al primo termine (bound superiore)."""
    m = len(periodogram)
    g = periodogram.max() / periodogram.sum()
    logp = np.log(m) + (m - 1) * np.log1p(-g)
    return float(g), float(min(1.0, np.exp(logp)))


def run(bits, label=""):
    b = bits[:CAP]
    N = len(b)
    x = b.astype(np.float64) - b.mean()
    if not x.any():
        return [], {"N": N, "nota": "sequenza costante"}
    L = int(min(MAX_LAG, N // 100))
    pvals, desc = [], {"N": N}
    if L >= 1:
        r = acf(x, L)
        lags = np.arange(1, L + 1)
        z = r * N / np.sqrt(N - lags)
        q = N * (N + 2) * np.sum(r ** 2 / (N - lags))
        p_lb = float(stats.chi2.sf(q, L))
        zmax = float(np.abs(z).max())
        p_one = 2 * stats.norm.sf(zmax)
        p_max = float(-np.expm1(L * np.log1p(-p_one))) if p_one < 1 else 1.0
        pvals += [
            {"test": "T4", "name": "ljung_box", "p": p_lb, "stat": float(q), "dof": L},
            {"test": "T4", "name": "acf_max_z", "p": p_max, "stat": zmax,
             "lag": int(lags[np.abs(z).argmax()])},
        ]
        desc.update(L=L, acf=r.tolist())
    spec = np.abs(np.fft.rfft(x)) ** 2 / N
    m = (N - 1) // 2
    per = spec[1:m + 1]
    g, p_g = fisher_g_pvalue(per)
    pvals.append({"test": "T4", "name": "fisher_g", "p": p_g, "stat": g,
                  "freq": float((per.argmax() + 1) / N)})
    # Spettro ridotto per i grafici: media su 1024 bande.
    bands = np.array_split(per / per.mean(), min(1024, len(per)))
    desc["spettro_bande"] = [float(c.mean()) for c in bands]
    desc["spettro_max_rel"] = float(per.max() / per.mean())
    return pvals, desc
