"""Caricamento dei dataset ed estrazione delle sequenze binarie.

Il formato interno dei file non è stato ispezionabile da questo ambiente (sorgenti bloccate
dalla rete). Per i dataset di Bell il caricatore usa quindi una mappatura esplicita delle
colonne, data/<nome>_mapping.json, da scrivere dopo `python -m qrand.pipeline inspect <nome>`.
Una mappatura sbagliata non passa inosservata: il sanity check (k = 196 per Delft, J > 0 e
p <= 1e-6 per NIST) blocca la pipeline.
"""

import json
import re
from pathlib import Path

import numpy as np

from . import DATA_DIR, RAW_DIR


class DatasetNonDisponibile(RuntimeError):
    pass


class FormatoNonVerificato(RuntimeError):
    pass


def status(name):
    path = RAW_DIR / name / "STATUS.json"
    if not path.exists():
        raise DatasetNonDisponibile(f"{name}: download non eseguito (make download)")
    st = json.loads(path.read_text())
    if st.get("esito") != "ok":
        raise DatasetNonDisponibile(f"{name}: {st.get('esito')} — {st.get('errore', '')}")
    return st


def data_files(name):
    return sorted(p for p in (RAW_DIR / name).rglob("*")
                  if p.is_file() and p.name not in ("STATUS.json", "links.txt"))


def inspect(name, max_lines=15):
    """Descrizione della struttura dei file grezzi, per scrivere la mappatura."""
    lines = []
    for p in data_files(name):
        lines.append(f"== {p.relative_to(RAW_DIR)} ({p.stat().st_size} byte)")
        if p.suffix.lower() in (".h5", ".hdf5"):
            import h5py
            with h5py.File(p, "r") as f:
                f.visititems(lambda k, v: lines.append(
                    f"   {k}: {getattr(v, 'shape', '')} {getattr(v, 'dtype', '')}"))
        else:
            try:
                with open(p, "r", errors="replace") as f:
                    for _, row in zip(range(max_lines), f):
                        lines.append("   " + row.rstrip()[:200])
            except OSError as e:
                lines.append(f"   (non leggibile: {e})")
    return "\n".join(lines)


def _mapping(name):
    path = DATA_DIR / f"{name}_mapping.json"
    if not path.exists():
        raise FormatoNonVerificato(
            f"{name}: manca {path.name}. Eseguire `python -m qrand.pipeline inspect {name}` e "
            f"scrivere la mappatura (vedi data/README.md).")
    return json.loads(path.read_text())


def _bell_from_table(name):
    """Mappatura tabellare: {"file": ..., "delimiter": ..., "skiprows": ..., "filter": {...},
    "columns": {"a_set": col, "b_set": col, "a_out": col, "b_out": col},
    "invert": [nomi di colonne da negare]}"""
    m = _mapping(name)
    path = RAW_DIR / name / m["file"]
    if path.suffix.lower() in (".h5", ".hdf5"):
        import h5py
        with h5py.File(path, "r") as f:
            cols = {k: np.asarray(f[v]).astype(np.int64).ravel() for k, v in m["columns"].items()}
    else:
        table = np.loadtxt(path, delimiter=m.get("delimiter"), skiprows=m.get("skiprows", 0),
                           dtype=float)
        cols = {k: table[:, int(v)].astype(np.int64) for k, v in m["columns"].items()}
        for col, value in m.get("filter", {}).items():
            keep = table[:, int(col)] == value
            cols = {k: v[keep] for k, v in cols.items()}
    for k in m.get("invert", []):
        cols[k] = 1 - cols[k]
    for k, v in cols.items():
        if not np.isin(v, (0, 1)).all():
            raise FormatoNonVerificato(f"{name}: la colonna {k} non è binaria")
    return cols


def bell_sequences(cols):
    a_set, b_set, a_out, b_out = (cols[k].astype(np.uint8)
                                  for k in ("a_set", "b_set", "a_out", "b_out"))
    seqs = {"A_set": a_set, "B_set": b_set, "A_out": a_out, "B_out": b_out,
            "coinc": a_out & b_out, "xor": a_out ^ b_out}
    joint = (a_set.astype(np.int64) << 3) | (b_set << 2) | (a_out << 1) | b_out
    return seqs, joint


def load_delft():
    status("delft")
    return _bell_from_table("delft")


def load_nist():
    status("nist")
    return _bell_from_table("nist")


def load_sycamore():
    """Serie per qubit e pesi di Hamming dal file di bitstring n53/m20 più grande."""
    status("sycamore")
    cands = [p for p in data_files("sycamore")
             if re.search(r"n53", p.name) and re.search(r"m20", p.name) and p.suffix == ".txt"]
    if not cands:
        raise FormatoNonVerificato(
            "sycamore: nessun file di testo con 'n53' e 'm20' nel nome; estrarre gli archivi e "
            "controllare con `inspect sycamore`")
    path = max(cands, key=lambda p: p.stat().st_size)
    rows = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not rows or any(len(r) != 53 or set(r) - {"0", "1"} for r in rows[:1000]):
        raise FormatoNonVerificato(f"sycamore: {path.name} non contiene righe di 53 caratteri 0/1")
    mat = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(len(rows), 53) - ord("0")
    seqs = {f"q{q:02d}": np.ascontiguousarray(mat[:, q]) for q in range(53)}
    return {"file": path.name, "n_bitstring": len(rows)}, seqs, mat.sum(axis=1)


def load_curby():
    """Bit dei round CURBy salvati; il campo con l'output casuale va indicato nella mappatura
    data/curby_mapping.json: {"field": "chiave.annidata", "encoding": "hex"|"base64"}."""
    import base64
    status("curby")
    m = _mapping("curby")
    out = []
    for p in data_files("curby"):
        node = json.loads(p.read_text())
        for key in m["field"].split("."):
            node = node[key]
        raw = bytes.fromhex(node) if m["encoding"] == "hex" else base64.b64decode(node)
        out.append(np.unpackbits(np.frombuffer(raw, dtype=np.uint8)))
    return np.concatenate(out) if out else np.zeros(0, np.uint8)
