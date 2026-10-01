# Pre-registrazione dell'analisi

**Stato:** fissata prima di eseguire la batteria su qualunque sequenza (dati reali o
controlli). Il commit git che introduce questo file precede ogni commit che contiene
risultati: la data del commit fa da marca temporale. Eventuali modifiche successive vanno
in fondo, nella sezione "Deviazioni", con data e motivo; il testo sopra non si modifica.

**Prima della stesura** sono stati eseguiti solo smoke test degli strumenti, su byte
casuali usa e getta che non fanno parte dell'analisi: compilazione e un'esecuzione
di `ea_non_iid` su 10^6 bit di `/dev/urandom` (H_min ≈ 0,84 bit/bit), una
misura dei tempi di `assess` (sts-2.1.2) su 10 × 10^6 bit e dei tempi di lzma e zstd.
Si è inoltre calcolato P(K ≥ 196 | Bin(245, 3/4)) = 0,0391 per fissare il criterio del
sanity check di Delft (sotto). Nessun test della batteria è stato eseguito.

## 1. Ipotesi

- **H0 (per ogni sequenza):** la sequenza è compatibile con il modello nullo indicato
  nella sezione 3. Il modello è i.i.d. tra prove successive, con la distribuzione
  marginale prevista dalla meccanica quantistica, i cui parametri si stimano dai dati.
- **H1:** esiste struttura sequenziale (dipendenza temporale, periodicità,
  ridondanza, "parole" ricorrenti) non spiegata dal modello nullo.

Le correlazioni *tra* Alice e Bob previste dalla MQ (la violazione di Bell stessa), il
bias dei rivelatori e la distribuzione di Porter–Thomas di Sycamore sono **struttura
attesa**. Non contano come anomalia; i modelli nulli sono costruiti per assorbirle.

## 2. Livelli e correzione per test multipli

- α = 0,01 per tutti i test.
- **Famiglia** = tutti i p-value confermativi di un dataset: tutte le sequenze e tutti
  i test. Ogni sequenza di controllo è una famiglia a sé, analizzata come se fosse un
  dataset.
- Correzione primaria: **Benjamini–Hochberg** con q = 0,01. Di confronto:
  **Bonferroni** con α = 0,01.
- Per i dati reali la **scoperta** si fa sulla prima metà di ogni sequenza (H1)
  con la correzione sulla famiglia di H1. Ogni p-value significativo dopo la correzione
  va poi **replicato** sulla seconda metà (H2), stesso test e stessi parametri, e conta
  come replicato se p_H2 < 0,01 e l'effetto ha lo stesso segno. Si chiama "anomalia"
  solo ciò che è significativo in H1 dopo la correzione *e* replicato in H2. L'analisi
  sulla lunghezza intera si riporta come descrittiva.
- I test con p-value discreto o approssimato restano tali e quali. Nessun p-value viene
  "aggiustato" a mano.

## 3. Sequenze e modelli nulli

| Sequenza | Definizione | Modello nullo |
|---|---|---|
| `A_set`, `B_set` | bit di scelta del setting per prova | i.i.d. Bernoulli(1/2) |
| `A_out`, `B_out` | esito per prova (Delft: lettura 0/1; NIST: 1 = rivelazione) | i.i.d. Bernoulli(p̂) |
| `coinc` | `A_out AND B_out` (coincidenza di esiti "1") | i.i.d. Bernoulli(p̂) |
| `xor` | `A_out XOR B_out` | i.i.d. Bernoulli(p̂) |
| `joint` | simbolo a 4 bit per prova (a_set, b_set, a_out, b_out) | prove i.i.d., distribuzione su 16 simboli libera |
| Sycamore `q00..q52` | serie temporale del qubit q nei campionamenti successivi di un file | i.i.d. Bernoulli(p̂_q) |
| Sycamore `hw` | peso di Hamming di ciascun bitstring | campionamenti i.i.d. |
| CURBy | bit di output pubblicati | i.i.d. Bernoulli(1/2) |
| Controlli | vedi sezione 6 | i.i.d. Bernoulli(1/2) |

`p̂` è la frequenza empirica degli "1" nella sequenza (o nella metà) analizzata.

## 4. Sanity check fisico (bloccante)

Si fa prima di qualunque altra analisi. Se fallisce, la pipeline si ferma per quel
dataset e il motivo va nel report.

- **Delft 2015 (Hensen et al., Nature 526, 682):** dai dati si ricalcolano il numero di
  prove n e il numero k di vittorie del gioco CHSH, poi il p-value
  P(K ≥ k | Bin(n, 3/4)) e S con errore standard. Il check **passa** se n = 245,
  k = 196, |p − 0,039| < 0,0005 e |S − 2,42| < 0,01. La convenzione dei segni del
  correlatore negativo è quella del paper; se non è ricavabile dalla documentazione si
  provano le 4 assegnazioni, si riportano tutte e si richiede che una sola soddisfi il
  criterio.
- **NIST 2015 (Shalm et al., PRL 115, 250402):** si ricalcolano i conteggi per
  combinazione di setting e la statistica CH/Eberhard J. Il check **passa** se J > 0 e
  se il p-value di un test (super)martingala con funzione di Bell fissata è
  ≤ 10^-6, cioè dello stesso ordine del valore dichiarato (p ≤ 2,3·10^-7 dopo
  l'aggiustamento per la predicibilità dei setting). La riproduzione esatta del valore
  PBR richiede il training set e la procedura del paper. Se non si riesce a stabilire
  la corrispondenza tra il file e la run del paper, il check è **"non riprodotto"** e
  la pipeline si ferma.
- **Sycamore, CURBy:** nessuna statistica di Bell. Si verificano solo il formato e il
  numero di bitstring o bit dichiarati dalla documentazione.

## 5. Batteria di test (lista chiusa)

Notazione: N = lunghezza in bit della sequenza analizzata.

### T1 — NIST SP 800-22 (suite C ufficiale sts-2.1.2, `assess`)
- Si applica solo alle sequenze con nullo Bernoulli(1/2) (setting, CURBy, controlli).
  Sulle altre è descrittiva: fallirebbe banalmente per il bias atteso.
- Segmentazione: n = 10^6 bit per sequenza, k = min(⌊N/10^6⌋, 100) sequenze dall'inizio
  (si usano quindi al più i primi 10^8 bit). Parametri di default di sts-2.1.2
  (M = 128, m = 9, m = 9, m = 10, m = 16, M = 500).
- p-value confermativi: se k ≥ 55, il P-value di uniformità di ogni riga di
  `finalAnalysisReport.txt` (188 righe; quelle con "----" si scartano). Se 1 ≤ k < 55,
  la combinazione di Fisher dei p-value di primo livello di ogni riga. Se N < 10^6,
  una sola sequenza di lunghezza N, solo con i test i cui requisiti minimi di
  lunghezza (SP 800-22 rev1a, §2) sono soddisfatti, usando i p-value di primo livello.
- Il test di proporzione NIST si riporta come descrittivo.

### T2 — NIST SP 800-90B (tool ufficiale `ea_non_iid`, modalità non-IID, 1 bit/simbolo)
- Si applica ai primi 10^6 simboli (requisito del tool). Con N < 10^6 non si applica.
- Riferimento: 19 sequenze surrogate i.i.d. Bernoulli(p̂) generate con `os.urandom`,
  stessa lunghezza, stesso tool.
- Esito **descrittivo** (è una stima, non un test): si riportano H_min(dati) e il rango
  tra i 19 riferimenti, p_rank = (1 + #{H_ref ≤ H_dati}) / 20. Si segnala la sequenza
  se H_min(dati) < min(H_ref). Il p_rank **non** entra nella famiglia, perché il suo
  minimo (0,05) è inutile dopo la correzione.

### T3 — Entropia a blocchi e test di ordine di Markov
- H_n per n = 1..20 su n-blocchi sovrapposti, con stimatore plug-in e correzione di
  Miller–Madow: Ĥ = Ĥ_plugin + (K̂ − 1)/(2 N_n ln 2), dove K̂ è il numero di blocchi con
  conteggio non nullo. h_n = H_{n+1} − H_n. Si mostrano in grafico con la banda di 5
  surrogati per permutazione. Valori con 2^n > N/10 sono marcati come non affidabili.
  Si usano i primi min(N, 2^27) bit.
- p-value confermativi: test del rapporto di verosimiglianza (G-test) tra catena di
  Markov di ordine k e modello i.i.d., con G = 2 Σ N(c,x) ln[N(c,x) N / (N(c) N(x))]
  ~ χ²(2^k − 1), per k ∈ {1, 2, 3, 4, 6, 8, 10, 12, 14, 16}, limitato ai k per cui
  N · min(p̂, 1−p̂)^(k+1) ≥ 5.

### T4 — Autocorrelazione e spettro
- Sequenza centrata x_i = b_i − p̂, primi min(N, 2^25) bit, lag 1..L con
  L = min(10^4, ⌊N/100⌋).
- p-value: (a) Ljung–Box Q_L ~ χ²(L); (b) max_k |z_k|, con z_k = r_k √(N−k) e
  correzione di Šidák su L lag.
- Spettro: periodogramma (FFT) della sequenza centrata. p-value: (c) test g di Fisher
  sul picco massimo (periodicità nascosta).

### T5 — Compressibilità
- Bit impacchettati in byte, primi min(N, 2^26) bit. Compressori: lzma
  (preset 9 | EXTREME, formato xz, nessun check) e zstd (livello 19).
- Nullo: M = 10 permutazioni casuali della stessa sequenza (stessa composizione).
- p-value (uno per compressore, unilaterale, "più comprimibile del nullo"):
  t = (C − media)/(sd·√(1 + 1/M)) ~ t(M−1).
- Descrittivo: C rispetto al limite di Shannon N·H(p̂)/8.

### T6 — Analisi "linguistica" (sul modello di Doyle & McCowan)
- **Parole fisse** di w ∈ {8, 16} bit, non sovrapposte, primi min(N, 2^27) bit.
  - (a) χ² di bontà di adattamento delle frequenze delle parole rispetto a
    Bernoulli(p̂)^w. Si fondono le celle con atteso < 5. Si fa solo se il numero di
    parole è ≥ 5 · 2^w.
  - (b) Pendenza di Zipf: regressione OLS di log f su log r per i ranghi 1..min(1000,
    tipi). Il p-value è bilaterale, con t di predizione rispetto a M = 10 surrogati per
    permutazione.
  - (c) Solo w = 8: G-test di indipendenza tra parole consecutive (tabella 256 × 256,
    entropia condizionale di ordine 1), se il numero di parole è ≥ 20 · 2^16.
  - Descrittivo: entropie condizionali di ordine 0..3 delle parole, con stimatore di
    Miller–Madow e banda dei surrogati.
- **Parole a lunghezza variabile:** token = run tra separatori. Separatore "1" (token =
  distanza tra "1" consecutivi) e, simmetricamente, separatore "0".
  - (d) χ² di adattamento della distribuzione delle lunghezze a Geometrica(p̂). Le
    lunghezze si troncano a L_max, scelta in modo che l'atteso sia ≥ 5.
  - (e) G-test di indipendenza tra lunghezze di token consecutivi (alfabeto troncato
    come sopra, celle con atteso medio ≥ 5).
  - Descrittivo: grafico rango-frequenza dei token.
- Per la sequenza `joint` (alfabeto di 16 simboli) e la sequenza `hw` di Sycamore si
  usa il G-test di ordine di Markov sui simboli (k = 1, 2, con lo stesso vincolo sugli
  attesi) al posto di T1, T3, T6, insieme a T4 e T5 sui simboli codificati.

## 6. Controlli (obbligatori)

Base comune: bit da `os.urandom`. Il modulo `secrets` usa la stessa sorgente del
kernel, quindi non si usa un secondo controllo separato.

- **Negativo forte:** `os.urandom`. N = 10^8 (1 replica, batteria completa con
  SP 800-22 100 × 10^6) e N = 10^7 (10 repliche, per il tasso empirico di falsi
  positivi).
- **Negativi deboli (PRNG difettosi):**
  - `lcg_low`: LCG con m = 2^32, a = 1664525, c = 1013904223, di cui si emette il byte
    meno significativo. I bit bassi hanno periodo corto, quindi ci si aspetta
    struttura evidente.
  - `randu_high`: RANDU (m = 2^31, a = 65539, c = 0), di cui si emette il byte più
    significativo. Il difetto è più sottile.
  - Lunghezze: come per il negativo forte (1 replica a 10^8, 10 repliche a 10^7).
- **Positivi:** testo ASCII (`data/controls/payload.txt`, testo della licenza
  Apache-2.0, 7 bit puri, ripetuto ciclicamente) inserito in bit da `os.urandom` a
  densità d ∈ {1%, 0,1%, 0,01%} dei bit, in due modalità:
  - `spread`: i bit del testo sostituiscono i bit nelle posizioni o + j·⌊1/d⌋
    (o è un offset casuale);
  - `block`: un blocco contiguo di ⌊d·N⌋ bit di testo in posizione casuale.
  - N = 10^8: 1 replica per configurazione, batteria completa. N = 10^7: 10 repliche
    per configurazione (T1 con 10 × 10^6, quindi combinazione di Fisher).
- **Potenza:**
  - Per test: frazione di repliche in cui almeno un p-value del test è < 0,01/K_test
    (Bonferroni dentro il test, dove K_test è il numero di p-value del test).
  - Per batteria: frazione di repliche con almeno una scoperta BH (q = 0,01) nella
    famiglia.
  - Si riporta la densità minima rilevata da ciascun test.

## 7. Regole di interpretazione (fissate)

- "Nessuna deviazione rilevata" vale solo insieme alla potenza dei controlli positivi
  alla stessa lunghezza N. Non equivale ad "assenza di struttura".
- Un p-value isolato, non replicato o non sopravvissuto alla correzione non si
  interpreta come segnale.
- Una deviazione replicata viene prima attribuita a cause note: deriva degli strumenti,
  bias del RNG, effetti di calibrazione o di lettura, struttura del formato dei file.
  Solo dopo si discutono altre ipotesi.

## 8. Sorgenti dei dati e politica di accesso

Si usano solo le pagine ufficiali: 4TU.ResearchData (Delft), il portale dati NIST
(Shalm), Dryad (Sycamore) e random.colorado.edu (CURBy). Non si usano copie di terze
parti. Se un dataset non è raggiungibile lo si segnala e si passa al successivo.

## Deviazioni

Tutte decise durante la validazione del codice su dati sintetici o di `os.urandom`, **prima**
di eseguire la batteria sui controlli o su dati reali.

1. **(2026-10-01) T4/T5 sulle sequenze di simboli (`joint`, Sycamore `hw`).** Il §5 prevedeva
   di applicarli "sui simboli codificati" in binario. Su prove sintetiche con statistica
   quantistica esatta questo dà p ≈ 0: i bit di una stessa prova sono correlati per
   costruzione (correlazioni di Bell, setting/esito), e centrare con la media globale crea
   una periodicità spuria di periodo pari alla larghezza del simbolo. Sotto il nullo "prove
   i.i.d." si usano quindi: T4 sui valori interi dei simboli centrati; T5 con un byte per
   simbolo e surrogati ottenuti permutando i simboli (non i bit).
2. **(2026-10-01) T1, dettagli d'implementazione di sts-2.1.2.**
   (a) Il P-value di uniformità si calcola in Python con la formula NIST (10 classi,
   igamc(9/2, χ²/2)), ma con atteso k/10 non troncato: il codice C usa una divisione intera.
   (b) La regola k ≥ 55 / Fisher si applica per riga al numero di sequenze *valide* della
   riga. Conta solo per Random Excursions (Variant), dove le sequenze con J < 500 non sono
   valide.
   (c) Se `assess` interrompe Random Excursions perché i cicli sono troppi
   (J > max(1000, n/100)), la riga viene esclusa e segnalata come "interrotta".
   (d) Con una sola sequenza (N < 2·10^6) ogni test gira in un processo `assess` separato.
   `assess` va in segfault nel rapporto finale quando un p-value non è in [0, 1], e questo
   può far perdere i risultati degli altri test. I p-value fuori da [0, 1] si scartano.
   (e) Il file d'ingresso è riempito fino a un multiplo di 4 byte, perché il lettore binario
   legge parole di 4 byte. Il riempimento non viene analizzato.
   (f) Per NonOverlappingTemplate, che non ha un minimo esplicito, si usa la soglia
   conservativa N ≥ 10^5.
   (g) `csrc.nist.gov` non è raggiungibile da questo ambiente. `tools/build_tools.sh` prova
   lo zip ufficiale e, se fallisce, usa la copia GitHub terrillmoore/NIST-Statistical-Test-Suite
   (codice NIST originale), con il commit fissato.
