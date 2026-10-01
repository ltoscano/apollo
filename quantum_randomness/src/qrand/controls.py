"""Sequenze di controllo: negativo forte (os.urandom), negativi deboli (LCG difettosi),
positivi (testo ASCII nascosto a bassa densità)."""

import os

import numpy as np

from . import DATA_DIR
from .bits import from_bytes, urandom_bits

PAYLOAD = DATA_DIR / "controls" / "payload.txt"
DENSITIES = (1e-2, 1e-3, 1e-4)
MODES = ("spread", "block")


def _rng():
    return np.random.default_rng(int.from_bytes(os.urandom(32), "little"))


def _lcg_bytes(n_bytes, m, a, c, shift, seed):
    """Byte (x >> shift) & 0xFF di un LCG x <- (a x + c) mod m, calcolato a blocchi con numpy."""
    out = np.empty(n_bytes, dtype=np.uint8)
    x = np.uint64(seed % m)
    # Salto di 2^16 passi: x_{i+j} = A_j x_i + C_j (mod m), tabelle per j = 1..block.
    block = 1 << 16
    A = np.empty(block, dtype=np.uint64)
    C = np.empty(block, dtype=np.uint64)
    a_j, c_j = 1, 0
    for j in range(block):
        a_j, c_j = (a * a_j) % m, (a * c_j + c) % m
        A[j], C[j] = a_j, c_j
    mask = np.uint64(m - 1)
    for start in range(0, n_bytes, block):
        k = min(block, n_bytes - start)
        xs = (A[:k] * x + C[:k]) & mask  # m potenza di 2: il wrap-around a 2^64 non altera il resto
        out[start:start + k] = ((xs >> np.uint64(shift)) & np.uint64(0xFF)).astype(np.uint8)
        x = xs[k - 1]
    return out


def lcg_low(n_bits, seed=None):
    """LCG m = 2^32, a = 1664525, c = 1013904223; byte meno significativo."""
    seed = int.from_bytes(os.urandom(4), "little") if seed is None else seed
    return from_bytes(_lcg_bytes((n_bits + 7) // 8, 2 ** 32, 1664525, 1013904223, 0, seed))[:n_bits]


def randu_high(n_bits, seed=None):
    """RANDU: m = 2^31, a = 65539, c = 0 (seme dispari); byte più significativo (bit 23..30)."""
    seed = (int.from_bytes(os.urandom(4), "little") | 1) if seed is None else seed
    return from_bytes(_lcg_bytes((n_bits + 7) // 8, 2 ** 31, 65539, 0, 23, seed))[:n_bits]


def payload_bits():
    return from_bytes(PAYLOAD.read_bytes())


def embed(cover, density, mode, payload=None):
    """Inserisce i bit del testo (ripetuto ciclicamente) nella sequenza di copertura."""
    bits = cover.copy()
    pay = payload_bits() if payload is None else payload
    N = len(bits)
    rng = _rng()
    if mode == "spread":
        stride = int(round(1 / density))
        offset = int(rng.integers(stride))
        pos = np.arange(offset, N, stride)
    elif mode == "block":
        length = int(density * N)
        start = int(rng.integers(0, N - length + 1))
        pos = np.arange(start, start + length)
    else:
        raise ValueError(mode)
    bits[pos] = np.resize(pay, pos.size)
    return bits, {"modo": mode, "densita": density, "bit_inseriti": int(pos.size),
                  "inizio": int(pos[0]) if pos.size else None}


def make(kind, n_bits, density=None, mode=None):
    if kind == "urandom":
        return urandom_bits(n_bits), {}
    if kind == "lcg_low":
        return lcg_low(n_bits), {}
    if kind == "randu_high":
        return randu_high(n_bits), {}
    if kind == "positivo":
        return embed(urandom_bits(n_bits), density, mode)
    raise ValueError(kind)
