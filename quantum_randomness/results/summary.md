# Riepilogo automatico dei risultati

Generato da `python -m qrand.pipeline report`. Non modificare a mano.

## Disponibilità dei dati

| Dataset | Fonte ufficiale | Esito download | Stato analisi | Licenza |
|---|---|---|---|---|
| delft | https://data.4tu.nl/articles/dataset/Loophole-free_Bell-inequality_violation_using_electron_spins_separated_by_1_3_kilometres/12703235 | non raggiungibile (URLError: <urlopen error Tunnel connection failed: 403 Forbidden>) | non analizzato | da leggere sulla pagina 4TU (non verificabile da questo ambiente) |
| nist | https://www.nist.gov/pml/applied-physics-division/bell-test-research-software-and-data/repository-bell-test-research-2 | non raggiungibile (URLError: <urlopen error Tunnel connection failed: 403 Forbidden>) | non analizzato | dati del governo federale USA (NIST); termini da leggere sulla pagina (non verificabile da questo ambiente) |
| sycamore | https://datadryad.org/dataset/doi:10.5061/dryad.k6t1rj8 | non raggiungibile (URLError: <urlopen error Tunnel connection failed: 403 Forbidden>) | non analizzato | CC0 1.0 (politica Dryad per tutti i dataset) |
| curby | https://random.colorado.edu/ | non raggiungibile (URLError: <urlopen error Tunnel connection failed: 403 Forbidden>) | non analizzato | non indicata nelle fonti consultate |

## Dataset

### delft

Non analizzato: delft: non raggiungibile — URLError: <urlopen error Tunnel connection failed: 403 Forbidden>

### nist

Non analizzato: nist: non raggiungibile — URLError: <urlopen error Tunnel connection failed: 403 Forbidden>

### sycamore

Non analizzato: sycamore: non raggiungibile — URLError: <urlopen error Tunnel connection failed: 403 Forbidden>

### curby

Non analizzato: curby: non raggiungibile — URLError: <urlopen error Tunnel connection failed: 403 Forbidden>


## Controlli a 10^8 bit (1 replica, batteria completa)

| Controllo | N bit | p-value nella famiglia | Scoperte BH (q=0,01) | Scoperte Bonferroni | Test con scoperte BH | H_min 90B (dati / min rif.) |
|---|---|---|---|---|---|---|

### p-value minimi per gruppo di test (grezzo, BH, Bonferroni sulla famiglia)

| Controllo | T1 p_min (p_BH, p_Bonf) | T3 p_min (p_BH, p_Bonf) | T4 p_min (p_BH, p_Bonf) | T5 p_min (p_BH, p_Bonf) | T6 p_min (p_BH, p_Bonf) |
|---|---|---|---|---|---|

## Tasso di falsi positivi (negativo forte)

- Repliche os.urandom a 10^7 bit: 4; p-value per famiglia: 211–211.
- Repliche con almeno una scoperta BH: 1 su 4 (scoperte totali: 1).
- Repliche segnalate per gruppo di test (Bonferroni nel gruppo): SP 800-22 1, SP 800-90B 0, Entropia/Markov 0, Autocorr./spettro 0, Compressione 0, Linguistica 0.

## Potenza (N = 10^7, 10 repliche per configurazione)

Rilevato = almeno un p-value del gruppo < 0,01/K (Bonferroni nel gruppo); per SP 800-90B, H_min sotto il minimo dei 19 riferimenti; batteria = almeno una scoperta BH nella famiglia.

| Configurazione (N = 10^7) | repliche | SP 800-22 | SP 800-90B | Entropia/Markov | Autocorr./spettro | Compressione | Linguistica | Batteria (BH) |
|---|---|---|---|---|---|---|---|---|
| os.urandom (negativo forte) | 4 | 25% | 0% | 0% | 0% | 0% | 0% | 25% |

### Densità minima rilevata con potenza ≥ 80%

| Modo | SP 800-22 | SP 800-90B | Entropia/Markov | Autocorr./spettro | Compressione | Linguistica | Batteria (BH) |
|---|---|---|---|---|---|---|---|
| spread | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato |
| block | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato | non rilevato |
