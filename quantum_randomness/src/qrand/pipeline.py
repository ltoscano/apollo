"""Punto d'ingresso della pipeline.

    python -m qrand.pipeline controls --set full|power [--workers N]
    python -m qrand.pipeline dataset delft|nist|sycamore|curby
    python -m qrand.pipeline inspect <dataset>
    python -m qrand.pipeline report
"""

import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from . import ALPHA, RESULTS_DIR, battery, bell, controls, loaders
from .bits import halves

CONTROL_SETS = {"full": {"n_bits": 10 ** 8, "replicates": 1},
                "power": {"n_bits": 10 ** 7, "replicates": 10}}
HEAVY_KEYS = ("acf", "spettro_bande", "rank_freq", "rank_freq_surrogato", "h_cond",
              "h_cond_surrogati", "h_surrogati", "H_ref", "stimatori")


def control_configs():
    cfg = [("urandom", {}), ("lcg_low", {}), ("randu_high", {})]
    for mode in controls.MODES:
        for d in controls.DENSITIES:
            cfg.append((f"pos_{mode}_{d:g}", {"density": d, "mode": mode}))
    return cfg


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def save(obj, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, default=_json_default, ensure_ascii=False))
    tmp.replace(path)


def slim(obj):
    if isinstance(obj, dict):
        return {k: slim(v) for k, v in obj.items() if k not in HEAVY_KEYS}
    if isinstance(obj, list):
        return [slim(v) for v in obj]
    return obj


def _control_job(set_name, name, params, rep, n_bits):
    t0 = time.time()
    if name in ("urandom", "lcg_low", "randu_high"):
        bits, info = controls.make(name, n_bits)
    else:
        bits, info = controls.make("positivo", n_bits, params["density"], params["mode"])
    res = battery.run_bits(bits, name, fair_null=True, log=lambda *a: None)
    fam = battery.correct([res])
    out = {"controllo": name, "set": set_name, "replica": rep, "parametri": params,
           "inserimento": info, "risultato": res,
           "famiglia": {k: v for k, v in fam.items() if k != "records"},
           "secondi": time.time() - t0}
    if rep > 0:
        out = slim(out)
    save(out, RESULTS_DIR / "controls" / set_name / f"{name}_r{rep}.json")
    return name, rep, fam["scoperte_bh"], time.time() - t0


def run_controls(set_name, workers):
    spec = CONTROL_SETS[set_name]
    jobs = []
    for name, params in control_configs():
        for rep in range(spec["replicates"]):
            if not (RESULTS_DIR / "controls" / set_name / f"{name}_r{rep}.json").exists():
                jobs.append((set_name, name, params, rep, spec["n_bits"]))
    print(f"controlli '{set_name}': {len(jobs)} job da eseguire", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_control_job, *j): j for j in jobs}
        for f in as_completed(futs):
            try:
                name, rep, nbh, sec = f.result()
                print(f"  {name} r{rep}: scoperte BH = {nbh} ({sec:.0f}s)", flush=True)
            except Exception:
                print(f"  ERRORE {futs[f][1:4]}:\n{traceback.format_exc()}", flush=True)


def _analyse(seqs, symbols, fair):
    """Batteria su H1 (scoperta), H2 (replica) e sequenza intera (descrittiva)."""
    parts = {}
    for part in ("H1", "H2", "intera"):
        results = []
        for label, bits in seqs.items():
            b = {"H1": halves(bits)[0], "H2": halves(bits)[1], "intera": bits}[part]
            results.append(battery.run_bits(b, label, fair_null=fair(label)))
        for label, (sym, alphabet) in symbols.items():
            s = {"H1": sym[: len(sym) // 2], "H2": sym[len(sym) // 2:], "intera": sym}[part]
            results.append(battery.run_symbols(s, alphabet, label))
        parts[part] = {"risultati": results, "famiglia": battery.correct(results)}
    # Replica: ogni scoperta BH in H1 si ripete in H2 (stesso test, stessa sequenza).
    h2 = {(r["sequenza"], r["test"], r["name"]): r for res in parts["H2"]["risultati"]
          for r in res["pvalues"]}
    anomalies = []
    for r in parts["H1"]["famiglia"]["records"]:
        if r["sig_bh"]:
            rep = h2.get((r["sequenza"], r["test"], r["name"]))
            anomalies.append({"sequenza": r["sequenza"], "test": r["test"], "nome": r["name"],
                              "p_H1": r["p"], "p_bh_H1": r["p_bh"],
                              "p_H2": rep["p"] if rep else None,
                              "replicata": bool(rep and rep["p"] < ALPHA)})
    return parts, anomalies


def run_dataset(name):
    out_dir = RESULTS_DIR / name
    try:
        if name in ("delft", "nist"):
            cols = loaders.load_delft() if name == "delft" else loaders.load_nist()
            check = (bell.delft_sanity if name == "delft" else bell.nist_sanity)(
                cols["a_set"], cols["b_set"], cols["a_out"], cols["b_out"])
            save(check, out_dir / "sanity.json")
            if not check["passato"]:
                print(f"{name}: sanity check NON superato — analisi interrotta (vedi sanity.json)")
                return 2
            seqs, joint = loaders.bell_sequences(cols)
            symbols = {"joint": (joint, 16)}
            fair = lambda label: label in ("A_set", "B_set")
        elif name == "sycamore":
            info, seqs, hw = loaders.load_sycamore()
            save({"passato": True, **info}, out_dir / "sanity.json")
            symbols = {"hw": (hw, 54)}
            fair = lambda label: False
        elif name == "curby":
            seqs, symbols = {"curby": loaders.load_curby()}, {}
            fair = lambda label: True
        else:
            raise ValueError(name)
    except (loaders.DatasetNonDisponibile, loaders.FormatoNonVerificato) as e:
        save({"esito": "non analizzato", "motivo": str(e)}, out_dir / "stato.json")
        print(f"{name}: non analizzato — {e}")
        return 1
    parts, anomalies = _analyse(seqs, symbols, fair)
    for part, content in parts.items():
        save(content, out_dir / f"analisi_{part}.json")
    save({"esito": "analizzato", "anomalie": anomalies}, out_dir / "stato.json")
    print(f"{name}: {len(anomalies)} scoperte in H1, "
          f"{sum(a['replicata'] for a in anomalies)} replicate in H2")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qrand.pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("controls")
    c.add_argument("--set", choices=CONTROL_SETS, required=True)
    c.add_argument("--workers", type=int, default=3)
    d = sub.add_parser("dataset")
    d.add_argument("name", choices=("delft", "nist", "sycamore", "curby"))
    i = sub.add_parser("inspect")
    i.add_argument("name", choices=("delft", "nist", "sycamore", "curby"))
    sub.add_parser("report")
    args = ap.parse_args(argv)
    if args.cmd == "controls":
        run_controls(args.set, args.workers)
    elif args.cmd == "dataset":
        return run_dataset(args.name)
    elif args.cmd == "inspect":
        print(loaders.inspect(args.name))
    elif args.cmd == "report":
        from . import report
        report.build()
    return 0


if __name__ == "__main__":
    sys.exit(main())
