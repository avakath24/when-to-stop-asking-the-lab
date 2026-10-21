#!/usr/bin/env bash
# Downloads the public Tox21 table (DeepChem copy) into ./data/tox21.csv
set -e
mkdir -p data
curl -L -o data/tox21.csv.gz https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/tox21.csv.gz
gunzip -f data/tox21.csv.gz
echo "Saved data/tox21.csv"
