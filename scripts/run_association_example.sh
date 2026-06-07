#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${ROOT}/.venv/bin/python"
if [[ ! -x "${PY}" ]]; then
  PY="python3"
fi

WEIGHT="${ROOT}/result/example_weight/weight_GRN.csv"
INFO="${ROOT}/result/example_info/info_GRN.csv"
GWAS="${ROOT}/example_data/gwas/AD_2019_tigar_GWAS.hg38.sorted.indexed.tsv.gz"
LD_PATTERN="${ROOT}/example_data/ld/CHR{chrom}_reference_cov.nochr.txt.gz"
GTF="${ROOT}/example_data/gene_example_55.bed"

if [[ ! -f "${WEIGHT}" || ! -f "${INFO}" ]]; then
  echo "Missing trained example weights. Run scripts/run_training_example.sh first." >&2
  exit 1
fi

if [[ ! -f "${GWAS}" || ! -f "${GWAS}.tbi" || ! -f "${ROOT}/example_data/ld/CHR1_reference_cov.nochr.txt.gz.tbi" ]]; then
  "${PY}" "${ROOT}/scripts/generate_example_association_data.py" \
    --weight "${WEIGHT}" \
    --gwas-dir "${ROOT}/example_data/gwas" \
    --ld-dir "${ROOT}/example_data/ld"
fi

mkdir -p "${ROOT}/result/example_association_ld" "${ROOT}/result/example_association_no_ld"

"${PY}" "${ROOT}/code/Association_GWAS/Association_GRNTWAS.py" \
  --gene_anno "${INFO}" \
  --weight "${WEIGHT}" \
  --Zscore "${GWAS}" \
  --LD_pattern "${LD_PATTERN}" \
  --window 100000 \
  --weight_threshold 0 \
  --thread 1 \
  --out_dir "${ROOT}/result/example_association_ld" \
  --out_twas_file "GRNTWAS_results_AD_2019_with_LD.tsv" \
  --gtf "${GTF}"

"${PY}" "${ROOT}/code/Association_GWAS/Association_GRNTWAS.py" \
  --gene_anno "${INFO}" \
  --weight "${WEIGHT}" \
  --Zscore "${GWAS}" \
  --window 100000 \
  --weight_threshold 0 \
  --thread 1 \
  --out_dir "${ROOT}/result/example_association_no_ld" \
  --out_twas_file "GRNTWAS_results_AD_2019_no_LD.tsv" \
  --gtf "${GTF}"

echo "Association examples complete:"
echo "  ${ROOT}/result/example_association_ld/GRNTWAS_results_AD_2019_with_LD.tsv"
echo "  ${ROOT}/result/example_association_no_ld/GRNTWAS_results_AD_2019_no_LD.tsv"
