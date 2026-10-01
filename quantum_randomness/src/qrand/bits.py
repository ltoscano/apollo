"""Utility per sequenze di bit (array numpy uint8 con valori 0/1)."""

import os

import numpy as np


def as_bits(x):
    b = np.asarray(x, dtype=np.uint8)
    if b.ndim != 1:
        raise ValueError("serve un array 1-D di bit")
    if b.size and b.max() > 1:
        raise ValueError("i bit devono valere 0 o 1")
    return b


def from_bytes(buf):
    return np.unpackbits(np.frombuffer(buf, dtype=np.uint8))


def to_bytes(bits):
    return np.packbits(as_bits(bits)).tobytes()


def urandom_bits(n):
    return from_bytes(os.urandom((n + 7) // 8))[:n]


def bernoulli_bits(n, p):
    """Bernoulli(p) i.i.d. con entropia da os.urandom (via numpy Generator seminato da urandom)."""
    if p == 0.5:
        return urandom_bits(n)
    rng = np.random.default_rng(int.from_bytes(os.urandom(32), "little"))
    return (rng.random(n) < p).astype(np.uint8)


def permuted(bits):
    rng = np.random.default_rng(int.from_bytes(os.urandom(32), "little"))
    return rng.permutation(bits)


def halves(bits):
    m = len(bits) // 2
    return bits[:m], bits[m:2 * m]


def block_codes(bits, n):
    """Codici interi (uint32) degli n-blocchi sovrapposti, n <= 32."""
    b = as_bits(bits)
    m = len(b) - n + 1
    if m <= 0:
        return np.zeros(0, dtype=np.uint32)
    codes = np.zeros(m, dtype=np.uint32)
    for i in range(n):
        codes <<= 1
        codes |= b[i:i + m]
    return codes


def iter_block_codes(bits, n_max):
    """Genera (n, codici degli n-blocchi sovrapposti) per n = 1..n_max, in modo incrementale."""
    b = as_bits(bits)
    codes = b.astype(np.uint32)
    yield 1, codes
    for n in range(2, n_max + 1):
        codes = (codes[:-1] << 1) | b[n - 1:]
        yield n, codes


def words(bits, w):
    """Parole non sovrapposte di w bit (w in {8, 16}) come interi, bit piu' significativo per primo."""
    b = as_bits(bits)
    packed = np.packbits(b[: (len(b) // w) * w])
    if w == 8:
        return packed.astype(np.int64)
    if w == 16:
        return packed.view(">u2").astype(np.int64)
    raise ValueError("w deve essere 8 o 16")


def binary_entropy(p):
    if p <= 0 or p >= 1:
        return 0.0
    return float(-p * np.log2(p) - (1 - p) * np.log2(1 - p))
