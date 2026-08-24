#!/usr/bin/env bash
set -euo pipefail

work="${1:-/tmp/sbatch-local-selftest}"
rm -rf "$work"
mkdir -p "$work/state" "$work/logs"

cat > "$work/normal.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=normal-test
#SBATCH --cpus-per-task=3
#SBATCH --mem=2g
#SBATCH --output=$work/logs/normal.%j.out
#SBATCH --error=$work/logs/normal.%j.err
set -euo pipefail
echo "job=\$SLURM_JOB_ID cpus=\$SLURM_CPUS_PER_TASK mem=\$SLURM_MEM_PER_NODE"
EOF

cat > "$work/array.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=array-test
#SBATCH --array=1-3
#SBATCH --cpus-per-task=2
#SBATCH --output=$work/logs/array.%A_%a.out
#SBATCH --error=$work/logs/array.%A_%a.err
set -euo pipefail
echo "array_job=\$SLURM_ARRAY_JOB_ID task=\$SLURM_ARRAY_TASK_ID"
EOF

VS_SBATCH_STATE_DIR="$work/state" sbatch "$work/normal.sh"
VS_SBATCH_STATE_DIR="$work/state" sbatch "$work/array.sh"

grep -R "cpus=3 mem=2g" "$work/logs" >/dev/null
grep -R "task=1" "$work/logs" >/dev/null
grep -R "task=2" "$work/logs" >/dev/null
grep -R "task=3" "$work/logs" >/dev/null

echo "PASS: local sbatch normal and array execution"
