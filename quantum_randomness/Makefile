# Pipeline completa: make all  (equivale a ./run_all.sh)
PY      ?= python3
WORKERS ?= 3
export PYTHONPATH := $(CURDIR)/src

.PHONY: all deps tools test download datasets controls controls-full controls-power report notebook clean-results

all: deps tools test download datasets controls report notebook

deps:
	$(PY) -m pip install -q -r requirements.txt

tools:
	./tools/build_tools.sh

test:
	$(PY) -m pytest -q tests

download:
	$(PY) data/download.py all

# Ogni dataset si ferma da solo (con stato registrato) se non è disponibile o se il
# sanity check fallisce; l'errore non interrompe il resto della pipeline.
datasets:
	-$(PY) -m qrand.pipeline dataset delft
	-$(PY) -m qrand.pipeline dataset nist
	-$(PY) -m qrand.pipeline dataset sycamore
	-$(PY) -m qrand.pipeline dataset curby

controls: controls-power controls-full

controls-full:
	$(PY) -m qrand.pipeline controls --set full --workers $(WORKERS)

controls-power:
	$(PY) -m qrand.pipeline controls --set power --workers $(WORKERS)

report:
	$(PY) -m qrand.pipeline report

notebook:
	$(PY) -m jupyter nbconvert --to notebook --execute notebooks/esplorazione.ipynb --output esplorazione.ipynb || true

clean-results:
	rm -rf results/controls results/delft results/nist results/sycamore results/curby results/figure results/summary.md
