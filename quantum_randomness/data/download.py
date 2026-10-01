#!/usr/bin/env python3
"""Scarica i dataset dalle sorgenti ufficiali elencate in sources.json.

Uso:  python data/download.py [delft|nist|sycamore|curby|all] [--nist-pattern REGEX]

Per ogni dataset scrive data/raw/<nome>/STATUS.json con esito, URL contattati, dimensioni e
SHA-256 dei file. Se una sorgente non è raggiungibile lo registra e passa alla successiva:
non si usano copie di terze parti.
"""

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
SOURCES = json.loads((HERE / "sources.json").read_text())
UA = {"User-Agent": "qrand-research/1.0 (riproducibilita' scientifica)"}


def fetch(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read(), r.headers.get("Content-Type", "")


def download(url, dest, timeout=120):
    dest.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 20):
            f.write(chunk)
            h.update(chunk)
    return {"file": dest.name, "url": url, "byte": dest.stat().st_size, "sha256": h.hexdigest()}


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def get_delft(args):
    files, _ = fetch(SOURCES["delft"]["api_file"])
    out = []
    for f in json.loads(files):
        out.append(download(f["download_url"], RAW / "delft" / f["name"]))
    return out


def get_nist(args):
    page = SOURCES["nist"]["pagina_ufficiale"]
    html, _ = fetch(page)
    p = Links()
    p.feed(html.decode("utf-8", errors="replace"))
    links = sorted({urllib.parse.urljoin(page, h) for h in p.links
                    if re.search(r"\.(h5|hdf5|zip|gz|bz2|dat|bin)(\?|$)", h, re.I)})
    (RAW / "nist").mkdir(parents=True, exist_ok=True)
    (RAW / "nist" / "links.txt").write_text("\n".join(links) + "\n")
    if not args.nist_pattern:
        print(f"NIST: {len(links)} link elencati in data/raw/nist/links.txt; "
              "scegliere la run del paper con --nist-pattern REGEX")
        return [{"elenco": links}]
    sel = [u for u in links if re.search(args.nist_pattern, u)]
    return [download(u, RAW / "nist" / Path(urllib.parse.urlparse(u).path).name) for u in sel]


def get_sycamore(args):
    api = SOURCES["sycamore"]["api"]
    meta = json.loads(fetch(api)[0])
    version = meta["_links"]["stash:version"]["href"]
    base = "https://datadryad.org"
    files = json.loads(fetch(base + version + "/files")[0])
    out = []
    for f in files["_embedded"]["stash:files"]:
        url = base + f["_links"]["stash:download"]["href"]
        out.append(download(url, RAW / "sycamore" / f["path"], timeout=600))
    return out


def get_curby(args):
    raw, ctype = fetch(SOURCES["curby"]["api"])
    dest = RAW / "curby" / f"round_latest_{int(time.time())}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(raw)
    return [{"file": dest.name, "byte": len(raw), "content_type": ctype,
             "sha256": hashlib.sha256(raw).hexdigest()}]


GETTERS = {"delft": get_delft, "nist": get_nist, "sycamore": get_sycamore, "curby": get_curby}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", nargs="?", default="all", choices=[*GETTERS, "all"])
    ap.add_argument("--nist-pattern", default=None)
    args = ap.parse_args()
    names = list(GETTERS) if args.dataset == "all" else [args.dataset]
    for name in names:
        status = {"dataset": name, "quando": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                  "fonte": SOURCES[name]["pagina_ufficiale"]}
        try:
            status["file"] = GETTERS[name](args)
            status["esito"] = "ok"
        except (urllib.error.URLError, OSError, KeyError, ValueError) as e:
            status["esito"] = "non raggiungibile"
            status["errore"] = f"{type(e).__name__}: {e}"
        (RAW / name).mkdir(parents=True, exist_ok=True)
        (RAW / name / "STATUS.json").write_text(json.dumps(status, indent=2, ensure_ascii=False))
        print(f"{name}: {status['esito']}" + (f" ({status.get('errore')})" if "errore" in status else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
