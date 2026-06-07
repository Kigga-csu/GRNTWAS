#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"${ROOT}/scripts/run_training_example.sh"
"${ROOT}/scripts/run_association_example.sh"

echo "Full GRNTWAS example pipeline complete."
