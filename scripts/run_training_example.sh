#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${ROOT}/.venv/bin/python"
if [[ ! -x "${PY}" ]]; then
  PY="python3"
fi

mkdir -p "${ROOT}/result/example_weight" "${ROOT}/result/example_info"

"${PY}" "${ROOT}/code/main.py" \
  --bed "${ROOT}/example_data/gene_example_55.bed" \
  --grn "${ROOT}/example_data/GRN_example.gexf" \
  --grn-tsv "${ROOT}/example_data/GRN_example.tsv" \
  --geno "${ROOT}/example_data/vcf/" \
  --exp "${ROOT}/example_data/tpm_normalized_gene_peer_vcf.csv" \
  --sample "${ROOT}/example_data/sample_ids.txt" \
  --eqtl "${ROOT}/example_data/TF_eQTL_example.txt" \
  --out-weight "${ROOT}/result/example_weight/" \
  --out-info "${ROOT}/result/example_info/" \
  --threads 1 \
  --method LTM \
  --mode lasso \
  --windows 100000 \
  --tf-numbers 10 \
  --maf 0.01 \
  --missing-rate 0.2 \
  --hwe 0.001 \
  --no-cv-r2 \
  --model elasticnet

echo "Training example complete:"
echo "  ${ROOT}/result/example_weight/weight_GRN.csv"
echo "  ${ROOT}/result/example_info/info_GRN.csv"
