import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qrand import autocorr, bell, controls, entropy, linguistic  # noqa: E402
from qrand.bits import block_codes, urandom_bits, words  # noqa: E402
from qrand.stats import benjamini_hochberg, bonferroni, fisher_combine  # noqa: E402


def test_bh_matches_reference():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    adj, rej = benjamini_hochberg(p, 0.05)
    # Calcolo a mano: min_{j>=i} p_(j) m / j.
    ref = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.074 * 10 / 7, 0.216, 0.216, 0.216]
    assert np.allclose(adj, ref, atol=1e-6)
    assert rej.sum() == 2
    assert bonferroni(p, 0.05)[1].sum() == 1


def test_fisher_uniform_is_not_significant():
    assert fisher_combine([0.5] * 10) > 0.5


def test_word_and_block_codes():
    b = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.uint8)
    assert words(b, 8).tolist() == [0b10110010, 0b11110000]
    assert block_codes(b, 3)[:3].tolist() == [0b101, 0b011, 0b110]


def test_markov_gtest_calibrated_under_null():
    ps = []
    for _ in range(200):
        _, tables = entropy.block_entropies(urandom_bits(4000), 2, markov_k=(2,))
        from qrand.stats import g_test_table
        ps.append(g_test_table(tables[2])[2])
    assert stats.kstest(ps, "uniform").pvalue > 0.001


def test_markov_detects_dependence():
    rng = np.random.default_rng(0)
    b = np.zeros(20000, dtype=np.uint8)
    b[0] = 1
    flips = rng.random(b.size) < 0.4  # catena di Markov: P(cambio) = 0.4
    for i in range(1, b.size):
        b[i] = b[i - 1] ^ flips[i]
    pv, _ = entropy.run(b)
    assert min(r["p"] for r in pv) < 1e-10


def test_autocorr_detects_period():
    b = urandom_bits(200_000)
    b[::50] = 1
    pv, _ = autocorr.run(b)
    assert {r["name"]: r["p"] for r in pv}["acf_max_z"] < 1e-6


def test_token_lengths_geometric():
    b = urandom_bits(10)
    b[:] = [1, 0, 0, 1, 1, 0, 1, 0, 0, 0]
    assert linguistic.token_lengths(b, 1).tolist() == [3, 1, 2]


def test_embed_density_and_payload():
    cover = np.zeros(1_000_000, dtype=np.uint8)
    for mode in controls.MODES:
        out, info = controls.embed(cover, 1e-3, mode)
        assert info["bit_inseriti"] == 1000
        assert out.sum() > 0
    pay = controls.payload_bits()
    assert pay.reshape(-1, 8)[:, 0].max() == 0  # ASCII a 7 bit


def test_chsh_reproduces_binomial_pvalue():
    assert stats.binom.sf(195, 245, 0.75) == pytest.approx(0.0391, abs=1e-4)
    rng = np.random.default_rng(3)
    x, y, a, b = bell.simulate_qm_trials(100_000, rng, visibility=1.0)
    conv, _, _ = bell.chsh_all_conventions(x, y, a, b)
    best = max(conv, key=lambda c: c["S"])
    assert best["meno_su"] == [1, 1]
    assert best["S"] == pytest.approx(2 * np.sqrt(2), abs=0.03)


def test_nist_sanity_rejects_local_data():
    rng = np.random.default_rng(4)
    n = 100_000
    x, y = rng.integers(0, 2, n), rng.integers(0, 2, n)
    a, b = rng.integers(0, 2, n), rng.integers(0, 2, n)
    assert not bell.nist_sanity(x, y, a, b)["passato"]
    x, y, a, b = bell.simulate_qm_trials(n, rng, visibility=0.9)
    assert bell.nist_sanity(x, y, a, b)["passato"]
