#!/usr/bin/env bash
# Scarica e compila gli strumenti NIST in tools/build/ (non versionata).
#   - SP 800-90B: repository ufficiale usnistgov/SP800-90B_EntropyAssessment (commit fissato).
#   - SP 800-22: suite C ufficiale sts-2.1.2. Si prova prima lo zip ufficiale NIST; se il sito
#     non è raggiungibile si usa la copia su GitHub terrillmoore/NIST-Statistical-Test-Suite,
#     che contiene il codice NIST originale (directory sts/) più correzioni minori di build.
# Dipendenze di sistema (Debian/Ubuntu): g++ make libbz2-dev libdivsufsort-dev libjsoncpp-dev
#   libssl-dev libgmp-dev libmpfr-dev
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p build && cd build

EA_REPO=https://github.com/usnistgov/SP800-90B_EntropyAssessment.git
EA_COMMIT=87c104d0ed4cbc96103e7b8b38d6f2c7e0a6b289
STS_ZIP=https://csrc.nist.gov/CSRC/media/Projects/Random-Bit-Generation/documents/sts-2_1_2.zip
STS_MIRROR=https://github.com/terrillmoore/NIST-Statistical-Test-Suite.git
STS_MIRROR_COMMIT=07cd57730402e9723b03f4775100773c36bfdc1b

if [ ! -x SP800-90B_EntropyAssessment/cpp/ea_non_iid ]; then
  [ -d SP800-90B_EntropyAssessment ] || git clone -q "$EA_REPO"
  git -C SP800-90B_EntropyAssessment fetch -q --depth 1 origin "$EA_COMMIT" 2>/dev/null || true
  git -C SP800-90B_EntropyAssessment checkout -q "$EA_COMMIT"
  make -C SP800-90B_EntropyAssessment/cpp non_iid
fi

if [ ! -x sts/assess ]; then
  rm -rf sts sts_src
  if curl -fsSL -m 60 -o sts.zip "$STS_ZIP" 2>/dev/null; then
    echo "sts-2.1.2: zip ufficiale NIST" > sts_provenance.txt
    unzip -q sts.zip -d sts_src
    src=$(dirname "$(find sts_src -name makefile -path '*sts*' | head -1)")
  else
    echo "sts-2.1.2: copia GitHub $STS_MIRROR @ $STS_MIRROR_COMMIT (csrc.nist.gov non raggiungibile)" > sts_provenance.txt
    git clone -q "$STS_MIRROR" sts_src
    git -C sts_src checkout -q "$STS_MIRROR_COMMIT"
    src=sts_src/sts
  fi
  mkdir -p "$src/obj"
  make -C "$src" >/dev/null
  mv "$src" sts
  rm -rf sts_src sts.zip
fi
echo "strumenti pronti:"; ls -la SP800-90B_EntropyAssessment/cpp/ea_non_iid sts/assess; cat sts_provenance.txt
