# Struttura non casuale in dati quantistici pubblici

Verifica riproducibile: i dati sperimentali quantistici pubblici (test di Bell di Delft e
NIST 2015, bitstring di Google Sycamore 2019, beacon CURBy) contengono struttura statistica
non prevista dalla meccanica quantistica?

- **Piano fissato in anticipo:** [`PREREGISTRATION.md`](PREREGISTRATION.md). Le deviazioni
  sono in fondo, con data e motivo.
- **Risultati e interpretazione:** [`REPORT.md`](REPORT.md).
- **Tabelle generate automaticamente:** `results/summary.md`, figure in `results/figure/`.

## Esecuzione

```bash
./run_all.sh          # oppure: make all
```

Il comando, in sequenza:

1. installa le dipendenze Python (`requirements.txt`);
2. compila gli strumenti NIST ufficiali in `tools/build/` (`tools/build_tools.sh`):
   `ea_non_iid` di SP 800-90B e `assess` di sts-2.1.2;
3. esegue i test di unità;
4. scarica i dati dalle sorgenti ufficiali (`data/download.py`);
5. per ogni dataset esegue il sanity check e poi la batteria (si ferma, registrando il
   motivo, se il dataset manca o se il check fallisce);
6. esegue i controlli (os.urandom, LCG difettosi, testo nascosto);
7. genera `results/summary.md`, le figure e il notebook.

Dipendenze di sistema per gli strumenti NIST (Debian/Ubuntu): `g++ make libbz2-dev
libdivsufsort-dev libjsoncpp-dev libssl-dev libgmp-dev libmpfr-dev`.

Il tempo è dominato dai controlli: circa 1–2 ore su 4 core, con circa 3 GB di RAM per
processo sulle sequenze da 10^8 bit (`make controls WORKERS=2` per ridurre la memoria).
I risultati già presenti vengono riutilizzati.

## Struttura

```
data/          download.py, sources.json (provenienza e licenze), README (mappature), controls/payload.txt
src/qrand/     bits, stats             utilità e statistica (BH, Bonferroni, G-test, Miller–Madow)
               nist22, sp90b           T1 SP 800-22 (sts-2.1.2), T2 SP 800-90B (ea_non_iid)
               entropy                 T3 H_n, h_n, G-test di ordine di Markov
               autocorr                T4 autocorrelazione e spettro
               compress                T5 lzma / zstd contro surrogati
               linguistic              T6 parole fisse e variabili, Zipf, entropie condizionali
               bell                    sanity check CHSH (Delft) e CH/Eberhard (NIST)
               loaders, controls       estrazione delle sequenze e generazione dei controlli
               battery, pipeline       esecuzione e correzione per test multipli; CLI
               report                  tabelle e figure
tests/         test di unità (pytest)
notebooks/     esplorazione.ipynb
results/       JSON dei risultati, summary.md, figure/
tools/         build_tools.sh
```
