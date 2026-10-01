"""T5 — Compressibilità (lzma, zstd) come proxy della complessità di Kolmogorov."""

import lzma

import numpy as np
import zstandard

from .bits import binary_entropy, permuted, to_bytes
from .stats import miller_madow_entropy, surrogate_t_pvalue

CAP = 2 ** 26
N_SURROGATES = 10


def size_lzma(buf):
    return len(lzma.compress(buf, format=lzma.FORMAT_XZ, check=lzma.CHECK_NONE,
                             preset=9 | lzma.PRESET_EXTREME))


def size_zstd(buf):
    return len(zstandard.ZstdCompressor(level=19).compress(buf))


COMPRESSORS = {"lzma": size_lzma, "zstd": size_zstd}


def run(bits, label="", symbols=False):
    """Con symbols=True l'input è una sequenza di simboli < 256: un byte per simbolo, surrogati
    per permutazione dei simboli, limite di Shannon dall'entropia empirica dei simboli."""
    if symbols:
        b = np.asarray(bits[: CAP // 8], dtype=np.uint8)
        N = len(b)
        data = b.tobytes()
        surrogates = [permuted(b).tobytes() for _ in range(N_SURROGATES)]
        shannon = N * miller_madow_entropy(np.bincount(b)) / 8
    else:
        b = bits[:CAP]
        N = len(b)
        data = to_bytes(b)
        surrogates = [to_bytes(permuted(b)) for _ in range(N_SURROGATES)]
        shannon = N * binary_entropy(float(b.mean())) / 8
    pvals, desc = [], {"N": N, "byte_grezzi": len(data), "limite_shannon": shannon}
    for name, fn in COMPRESSORS.items():
        c = fn(data)
        cs = [fn(s) for s in surrogates]
        t, pv = surrogate_t_pvalue(c, cs, "less")
        pvals.append({"test": "T5", "name": name, "p": pv, "stat": t})
        desc[name] = {"dimensione": c, "surrogati": cs, "rapporto": c / len(data),
                      "rapporto_shannon": c / shannon if shannon else None}
    return pvals, desc
