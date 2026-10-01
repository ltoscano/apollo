"""Tabelle e figure dai risultati in results/: results/summary.md e results/figure/*.png."""

import json
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import ALPHA, DATA_DIR, RAW_DIR, RESULTS_DIR  # noqa: E402
from .controls import DENSITIES, MODES  # noqa: E402

FIG = RESULTS_DIR / "figure"
TESTS = ["T1", "T2", "T3", "T4", "T5", "T6"]
TEST_NAMES = {"T1": "SP 800-22", "T2": "SP 800-90B", "T3": "Entropia/Markov",
              "T4": "Autocorr./spettro", "T5": "Compressione", "T6": "Linguistica"}
# Palette categoriale di riferimento (ordine fisso), superficie e inchiostri.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
TEST_COLOR = {"T1": SERIES[0], "T3": SERIES[1], "T4": SERIES[2], "T5": SERIES[3],
              "T6": SERIES[4], "batteria": SERIES[6]}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "font.size": 10, "legend.frameon": False,
    "savefig.dpi": 130, "savefig.bbox": "tight",
})

CONTROL_LABEL = {"urandom": "os.urandom (negativo forte)", "lcg_low": "LCG byte basso (negativo debole)",
                 "randu_high": "RANDU byte alto (negativo debole)"}


def label(name):
    if name in CONTROL_LABEL:
        return CONTROL_LABEL[name]
    _, mode, d = name.split("_")
    return f"testo ASCII, {mode}, d = {float(d) * 100:g}%"


def load_controls(set_name):
    out = defaultdict(list)
    for p in sorted((RESULTS_DIR / "controls" / set_name).glob("*.json")):
        r = json.loads(p.read_text())
        out[r["controllo"]].append(r)
    return out


def detected_by_test(run):
    """{test: bool} con Bonferroni dentro il test (p < 0.01 / K_test) e "batteria" (BH)."""
    recs = [r for r in run["risultato"]["pvalues"] if r.get("confermativo")]
    det = {}
    for t in TESTS:
        ps = [r["p"] for r in recs if r["test"] == t and r["p"] == r["p"]]
        det[t] = bool(ps) and min(ps) < ALPHA / len(ps)
    t2 = run["risultato"]["descrittivi"].get("T2", {})
    det["T2"] = bool(t2.get("segnalata"))
    det["batteria"] = run["famiglia"]["scoperte_bh"] > 0
    return det


def config_order():
    names = ["urandom", "lcg_low", "randu_high"]
    names += [f"pos_{m}_{d:g}" for m in MODES for d in DENSITIES]
    return names


def availability_table():
    sources = json.loads((DATA_DIR / "sources.json").read_text())
    rows = ["| Dataset | Fonte ufficiale | Esito download | Stato analisi | Licenza |",
            "|---|---|---|---|---|"]
    for name, src in sources.items():
        st_path = RAW_DIR / name / "STATUS.json"
        st = json.loads(st_path.read_text()) if st_path.exists() else {"esito": "non eseguito"}
        an_path = RESULTS_DIR / name / "stato.json"
        an = json.loads(an_path.read_text()) if an_path.exists() else {"esito": "—"}
        err = st.get("errore", "")
        rows.append(f"| {name} | {src['pagina_ufficiale']} | {st['esito']}"
                    f"{' (' + err + ')' if err else ''} | {an['esito']} | {src['licenza']} |")
    return "\n".join(rows)


def full_controls_table(full):
    rows = ["| Controllo | N bit | p-value nella famiglia | Scoperte BH (q=0,01) | Scoperte Bonferroni "
            "| Test con scoperte BH | H_min 90B (dati / min rif.) |",
            "|---|---|---|---|---|---|---|"]
    for name in config_order():
        if name not in full:
            continue
        r = full[name][0]
        recs = [x for x in r["risultato"]["pvalues"] if x.get("confermativo")]
        p = np.array([x["p"] for x in recs])
        from .stats import benjamini_hochberg, bonferroni
        _, bh = benjamini_hochberg(p, ALPHA)
        _, bf = bonferroni(p, ALPHA)
        tests = sorted({recs[i]["test"] for i in np.flatnonzero(bh)})
        t2 = r["risultato"]["descrittivi"].get("T2", {})
        h = f"{t2['H_min']:.3f} / {t2['H_ref_min']:.3f}" if t2.get("applicabile") else "—"
        rows.append(f"| {label(name)} | {r['risultato']['N']:.0e} | {len(p)} | {int(bh.sum())} | "
                    f"{int(bf.sum())} | {', '.join(tests) or '—'} | {h} |")
    return "\n".join(rows)


def pvalue_table(full):
    """Per ogni controllo a 10^8 bit e ogni gruppo di test: p minimo grezzo e corretto."""
    rows = ["| Controllo | " + " | ".join(f"{t} p_min (p_BH, p_Bonf)" for t in TESTS if t != "T2") + " |",
            "|---|" + "---|" * 5]
    from .stats import benjamini_hochberg, bonferroni
    for name in config_order():
        if name not in full:
            continue
        recs = [x for x in full[name][0]["risultato"]["pvalues"] if x.get("confermativo")]
        p = np.array([x["p"] for x in recs])
        bh, _ = benjamini_hochberg(p, ALPHA)
        bf, _ = bonferroni(p, ALPHA)
        cells = []
        for t in TESTS:
            if t == "T2":
                continue
            idx = [i for i, x in enumerate(recs) if x["test"] == t]
            if not idx:
                cells.append("—")
                continue
            i = min(idx, key=lambda j: p[j])
            cells.append(f"{p[i]:.2g} ({bh[i]:.2g}, {bf[i]:.2g})")
        rows.append(f"| {label(name)} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def power_table(power):
    keys = TESTS + ["batteria"]
    rows = ["| Configurazione (N = 10^7) | repliche | " + " | ".join(
        TEST_NAMES.get(k, "Batteria (BH)") for k in keys) + " |", "|---|---|" + "---|" * len(keys)]
    rates = {}
    for name in config_order():
        runs = power.get(name, [])
        if not runs:
            continue
        dets = [detected_by_test(r) for r in runs]
        rates[name] = {k: float(np.mean([d[k] for d in dets])) for k in keys}
        rows.append(f"| {label(name)} | {len(runs)} | " + " | ".join(
            f"{rates[name][k]:.0%}" for k in keys) + " |")
    return "\n".join(rows), rates


def min_density_table(rates):
    keys = TESTS + ["batteria"]
    rows = ["| Modo | " + " | ".join(TEST_NAMES.get(k, "Batteria (BH)") for k in keys) + " |",
            "|---|" + "---|" * len(keys)]
    for mode in MODES:
        cells = []
        for k in keys:
            ok = [d for d in DENSITIES if rates.get(f"pos_{mode}_{d:g}", {}).get(k, 0) >= 0.8]
            cells.append(f"{min(ok) * 100:g}%" if ok else "non rilevato")
        rows.append(f"| {mode} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def fig_power(rates):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, mode in zip(axes, MODES):
        x = [d * 100 for d in DENSITIES]
        for k in ["T1", "T3", "T4", "T5", "T6", "batteria"]:
            y = [rates.get(f"pos_{mode}_{d:g}", {}).get(k, np.nan) for d in DENSITIES]
            ax.plot(x, y, marker="o", markersize=6, color=TEST_COLOR[k],
                    label=TEST_NAMES.get(k, "Batteria (BH)"), linestyle="--" if k == "batteria" else "-")
        ax.set_xscale("log")
        ax.set_xlabel("densità del testo inserito (% dei bit)")
        ax.set_title(f"modo '{mode}'", color=INK)
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("frazione di repliche rilevate")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle("Potenza di rilevazione dei controlli positivi (N = 10^7 bit, 10 repliche)")
    fig.savefig(FIG / "potenza.png")
    plt.close(fig)


def fig_hn(full):
    names = [n for n in ["urandom", "randu_high", "lcg_low", "pos_block_0.01", "pos_spread_0.01"] if n in full]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for i, name in enumerate(names):
        d = full[name][0]["risultato"]["descrittivi"]["T3"]
        h = np.array(d["h"])
        n = np.arange(len(h))
        ok = np.array(d["affidabile"])
        ax.plot(n[ok], h[ok], marker="o", markersize=4, color=SERIES[i], label=label(name))
        if (~ok).any():
            ax.plot(n[~ok], h[~ok], linestyle=":", color=SERIES[i])
        if name == "urandom" and d.get("h_surrogati"):
            band = np.array(d["h_surrogati"])
            ax.fill_between(n, band.min(0), band.max(0), color=GRID, alpha=0.8, label="surrogati (permutazione)")
    ax.set_xlabel("n (lunghezza del contesto)")
    ax.set_ylabel("h_n = H_{n+1} − H_n  [bit]")
    ax.set_title("Entropia condizionale h_n (Miller–Madow), controlli a 10^8 bit")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "hn_controlli.png")
    # dettaglio vicino a 1 bit
    ax.set_ylim(0.9990, 1.0003)
    ax.set_title("h_n — dettaglio vicino a 1 bit")
    fig.savefig(FIG / "hn_controlli_zoom.png")
    plt.close(fig)


def fig_zipf(full):
    names = [n for n in ["urandom", "randu_high", "lcg_low", "pos_block_0.01"] if n in full]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, w in zip(axes, (8, 16)):
        for i, name in enumerate(names):
            d = full[name][0]["risultato"]["descrittivi"]["T6"][f"parole{w}"]
            rf = np.array(d["rank_freq"], dtype=float)
            ax.plot(np.arange(1, rf.size + 1), rf / rf.sum(), color=SERIES[i], label=label(name))
        d = full["urandom"][0]["risultato"]["descrittivi"]["T6"][f"parole{w}"]
        rf = np.array(d["rank_freq_surrogato"], dtype=float)
        ax.plot(np.arange(1, rf.size + 1), rf / rf.sum(), color=INK2, linestyle=":",
                label="surrogato urandom (permutazione)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("rango")
        ax.set_ylabel("frequenza relativa (primi 1000 ranghi)")
        ax.set_title(f"Rango-frequenza, parole di {w} bit")
    axes[1].legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.savefig(FIG / "zipf_controlli.png")
    plt.close(fig)


def fig_acf(full):
    names = [n for n in ["urandom", "randu_high", "pos_spread_0.01"] if n in full]
    fig, axes = plt.subplots(len(names), 1, figsize=(9, 2.2 * len(names)), sharex=True)
    for ax, (i, name) in zip(np.atleast_1d(axes), enumerate(names)):
        d = full[name][0]["risultato"]["descrittivi"]["T4"]
        r = np.array(d["acf"])
        N = d["N"]
        ax.plot(np.arange(1, r.size + 1), r * np.sqrt(N), color=SERIES[i], linewidth=0.8)
        ax.axhline(4.5, color=INK2, linestyle=":", linewidth=1)
        ax.axhline(-4.5, color=INK2, linestyle=":", linewidth=1)
        ax.set_ylabel("z ≈ r_k √N")
        ax.set_title(label(name), fontsize=9)
    np.atleast_1d(axes)[-1].set_xlabel("lag")
    fig.suptitle("Autocorrelazione standardizzata (lag 1..10^4); linee: ±4,5 (≈ soglia Šidák)")
    fig.savefig(FIG / "acf_controlli.png")
    plt.close(fig)


def fig_pvalue_hist(full):
    if "urandom" not in full:
        return
    recs = [x for x in full["urandom"][0]["risultato"]["pvalues"] if x.get("confermativo")]
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.hist([x["p"] for x in recs], bins=20, range=(0, 1), color=SERIES[0], edgecolor=SURFACE)
    ax.axhline(len(recs) / 20, color=INK2, linestyle=":", label="atteso sotto H0")
    ax.set_xlabel("p-value")
    ax.set_ylabel("conteggio")
    ax.set_title(f"Calibrazione: {len(recs)} p-value della famiglia os.urandom (10^8 bit)")
    ax.legend()
    fig.savefig(FIG / "calibrazione_pvalue.png")
    plt.close(fig)


def type1_table(power):
    runs = power.get("urandom", [])
    if not runs:
        return "(nessuna replica)"
    dets = [detected_by_test(r) for r in runs]
    n_disc = [r["famiglia"]["scoperte_bh"] for r in runs]
    m = [r["famiglia"]["m"] for r in runs]
    lines = [f"- Repliche os.urandom a 10^7 bit: {len(runs)}; p-value per famiglia: {min(m)}–{max(m)}.",
             f"- Repliche con almeno una scoperta BH: {sum(x > 0 for x in n_disc)} su {len(runs)} "
             f"(scoperte totali: {sum(n_disc)}).",
             "- Repliche segnalate per gruppo di test (Bonferroni nel gruppo): " + ", ".join(
                 f"{TEST_NAMES[t]} {sum(d[t] for d in dets)}" for t in TESTS) + "."]
    return "\n".join(lines)


def dataset_sections():
    out = []
    for name in ("delft", "nist", "sycamore", "curby"):
        st_path = RESULTS_DIR / name / "stato.json"
        if not st_path.exists():
            continue
        st = json.loads(st_path.read_text())
        out.append(f"### {name}\n")
        if st["esito"] != "analizzato":
            out.append(f"Non analizzato: {st.get('motivo', '')}\n")
            continue
        san = json.loads((RESULTS_DIR / name / "sanity.json").read_text())
        out.append(f"Sanity check: {'superato' if san.get('passato') else 'NON superato'}.\n")
        h1 = json.loads((RESULTS_DIR / name / "analisi_H1.json").read_text())["famiglia"]
        out.append(f"H1: {h1['m']} p-value, scoperte BH {h1['scoperte_bh']}, Bonferroni {h1['scoperte_bonf']}.\n")
        rows = ["| Sequenza | Test | p (H1) | p_BH | p_Bonf |", "|---|---|---|---|---|"]
        for r in sorted(h1["records"], key=lambda x: x["p"])[:15]:
            rows.append(f"| {r['sequenza']} | {r['test']}:{r['name']} | {r['p']:.3g} | {r['p_bh']:.3g} | {r['p_bonf']:.3g} |")
        out.append("\n".join(rows) + "\n")
        for a in st["anomalie"]:
            out.append(f"- Scoperta {a['sequenza']} {a['test']}:{a['nome']}: p_H1 = {a['p_H1']:.3g}, "
                       f"p_H2 = {a['p_H2']}, replicata = {a['replicata']}\n")
    return "\n".join(out)


def build():
    FIG.mkdir(parents=True, exist_ok=True)
    full, power = load_controls("full"), load_controls("power")
    ptab, rates = power_table(power)
    for f in (fig_hn, fig_zipf, fig_acf, fig_pvalue_hist):
        if full:
            f(full)
    if rates:
        fig_power(rates)
    md = ["# Riepilogo automatico dei risultati", "",
          "Generato da `python -m qrand.pipeline report`. Non modificare a mano.", "",
          "## Disponibilità dei dati", "", availability_table(), "",
          "## Dataset", "", dataset_sections(), "",
          "## Controlli a 10^8 bit (1 replica, batteria completa)", "", full_controls_table(full), "",
          "### p-value minimi per gruppo di test (grezzo, BH, Bonferroni sulla famiglia)", "",
          pvalue_table(full), "",
          "## Tasso di falsi positivi (negativo forte)", "", type1_table(power), "",
          "## Potenza (N = 10^7, 10 repliche per configurazione)", "",
          "Rilevato = almeno un p-value del gruppo < 0,01/K (Bonferroni nel gruppo); per SP 800-90B, "
          "H_min sotto il minimo dei 19 riferimenti; batteria = almeno una scoperta BH nella famiglia.", "",
          ptab, "", "### Densità minima rilevata con potenza ≥ 80%", "", min_density_table(rates), ""]
    (RESULTS_DIR / "summary.md").write_text("\n".join(md))
    print("scritto results/summary.md e results/figure/")
