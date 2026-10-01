# Dati

I dati **non** sono versionati: `data/download.py` li scarica in `data/raw/<dataset>/`,
che è esclusa da git. Per ogni dataset lo script scrive `data/raw/<dataset>/STATUS.json`
con l'esito del download, gli URL contattati, la dimensione e lo SHA-256 di ogni file.
Provenienza, formato atteso e licenza sono in `sources.json`.

```bash
python data/download.py all                     # tutti i dataset
python data/download.py nist                    # elenca i link della pagina NIST in raw/nist/links.txt
python data/download.py nist --nist-pattern '...regex della run del paper...'
```

## Sorgenti (solo pagine ufficiali)

| Dataset | Pagina | Contenuto atteso | Licenza |
|---|---|---|---|
| Delft 2015 (Hensen et al.) | 4TU.ResearchData, articolo 12703235 | per ogni prova: setting (RNG) e lettura in A e B; n = 245 | da leggere sulla pagina |
| NIST 2015 (Shalm et al.) | nist.gov, *Bell Test Data – processed compressed* | HDF5 per run (settembre 2015) | dati del governo USA; da leggere sulla pagina |
| Sycamore 2019 (Arute et al.) | Dryad doi:10.5061/dryad.k6t1rj8 | bitstring osservati per circuito (2,36 GB) | CC0 (politica Dryad) |
| CURBy | random.colorado.edu | JSON dell'ultimo round del beacon | non indicata |

Nell'ambiente in cui è stato preparato il repository, tutte e quattro le sorgenti
rispondono **403 dal proxy di rete** (`STATUS.json`: "non raggiungibile"). Formato e
numero di eventi non sono quindi verificati sui file reali: quanto sopra viene dagli
articoli e dalle pagine indicizzate.

## Mappature delle colonne (Delft, NIST, CURBy)

Non avendo potuto vedere i file, i caricatori **non indovinano** la struttura. Dopo il
download:

```bash
PYTHONPATH=src python -m qrand.pipeline inspect delft     # elenca file, intestazioni, dataset HDF5
```

poi si scrive `data/delft_mapping.json` (analogo `nist_mapping.json`):

```json
{
  "file": "nome_del_file.txt",
  "delimiter": null,
  "skiprows": 1,
  "filter": {"7": 1},
  "columns": {"a_set": 2, "b_set": 3, "a_out": 4, "b_out": 5},
  "invert": []
}
```

- `columns`: indici di colonna (file di testo) o percorsi dei dataset (HDF5) per setting
  ed esito di Alice e Bob.
- `filter`: colonna → valore, per tenere solo le prove event-ready usate nel paper.
- `invert`: colonne da negare, se la convenzione 0/1 è opposta.

Per CURBy, `curby_mapping.json` è `{"field": "chiave.annidata", "encoding": "hex"}`.

Una mappatura sbagliata non passa inosservata. Il sanity check deve riprodurre n = 245,
k = 196, p = 0,039 e S = 2,42 (Delft), oppure J > 0 e p ≤ 10^-6 (NIST); se non ci riesce,
la pipeline si ferma.

## Controlli

`controls/payload.txt` è il testo della licenza Apache-2.0 (ASCII a 7 bit, 11 425 byte),
usato come messaggio nascosto nei controlli positivi.
