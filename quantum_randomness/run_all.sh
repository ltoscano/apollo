#!/usr/bin/env bash
# Esegue l'intera pipeline: dipendenze, strumenti NIST, test di unità, download, analisi dei
# dataset, controlli, report. I risultati già calcolati vengono riutilizzati (make è
# idempotente per i controlli); `make clean-results` per ripartire da zero.
set -euo pipefail
cd "$(dirname "$0")"
make all "$@"
