#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(cd "$script_dir/.." && pwd)/config/run.env"

required_paths=(
    "$VS_DATABASE_DIR/ref/GRCh38.p14"
    "$VS_DATABASE_DIR/nr"
    "$VS_DATABASE_DIR/nt"
    "$VS_DATABASE_DIR/core_nt"
    "$VS_DATABASE_DIR/nr_dmnd/nr.dmnd"
    "$VS_DATABASE_DIR/core_nt_mmseqs/core_nt.dbtype"
    "$VS_DATABASE_DIR/VirusDBNR/virus_nr.dmnd"
    "$VS_DATABASE_DIR/VirusDBNT/virus_nt.dbtype"
    "$VS_DATABASE_DIR/taxdump"
    "$VS_DATABASE_DIR/famdb"
    "$VS_DATABASE_DIR/vhunter_acc.db"
    "$VS_DATABASE_DIR/taxa.sqlite"
    "$VS_CONFIG_FILE"
    "$VS_CFG_FILE"
)

for path in "${required_paths[@]}"; do
    [[ -e "$path" ]] || {
        echo "ERROR: missing $path" >&2
        exit 1
    }
done

if [[ "$VS_DB_PROFILE" == "vs" ]]; then
    compgen -G "$VS_DATABASE_DIR/famdb/*.h5" >/dev/null || {
        echo "ERROR: FamDB files are missing." >&2
        exit 1
    }
else
    [[ -s "$VS_DATABASE_DIR/test/repeat_library.fa" ]] || {
        echo "ERROR: test RepeatMasker library is missing." >&2
        exit 1
    }
fi

grep -q '^vhunter = ' "$VS_CFG_FILE"
grep -q '^ncbi_taxadb = ' "$VS_CFG_FILE"

if [[ "$VS_DB_PROFILE" == "test" ]]; then
    [[ -z "${VS_EXTRA_FLAGS:-}" ]]
fi

case "$VS_EXECUTION_MODE" in
    docker)
        docker image inspect "$VS_IMAGE" >/dev/null
        docker run --rm -i \
            --platform linux/amd64 \
            --entrypoint /bin/bash \
            -e FAMDB_DATA_DIR=/database/famdb \
            -v "$VS_DATABASE_DIR:/database:ro" \
            --mount "type=bind,src=$VS_CONFIG_FILE,dst=/opt/virusseeker/scripts/VS.config,readonly" \
            --mount "type=bind,src=$VS_CONFIG_FILE,dst=/opt/virusseeker/scripts/VS_Read_Counter/VS.config,readonly" \
            --mount "type=bind,src=$VS_CFG_FILE,dst=/opt/virusseeker/scripts/VS.cfg,readonly" \
            "$VS_IMAGE" \
            -lc '
                set -Eeuo pipefail
                /opt/conda/envs/vs/bin/python -c "from Bio import SeqIO; print(\"Biopython SeqIO OK\")"
                /opt/conda/envs/vs/bin/python - <<PY
import configparser
cfg = configparser.ConfigParser()
cfg.read("/opt/virusseeker/scripts/VS.cfg")
assert cfg["paths"]["vhunter"] == "/database/vhunter_acc.db"
assert cfg["paths"]["ncbi_taxadb"] == "/database/taxa.sqlite"
print("VS.cfg Docker paths OK")
PY
                blastdbcmd -info -db /database/nt/nt >/dev/null
                sbatch --version
            '
        ;;

    slurm)
        command -v sbatch >/dev/null
        sbatch --version
        [[ -r "$VS_NATIVE_VS_HOME/queue.sh" ]]
        bash -n "$VS_NATIVE_VS_HOME/queue.sh"
        ;;

    *)
        echo "ERROR: invalid VS_EXECUTION_MODE=$VS_EXECUTION_MODE" >&2
        exit 1
        ;;
esac

echo "PASS: VirusSeeker $VS_EXECUTION_MODE installation is valid."
