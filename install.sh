#!/usr/bin/env bash
set -Eeuo pipefail

# VirusSeeker 2.0 dual-runtime installer.
#
# Runtime modes:
#   docker  - VirusSeeker runs in the Docker image and uses the local sbatch shim.
#   slurm   - VirusSeeker runs natively from the Git checkout and uses the site's
#             real sbatch/SLURM installation.
#
# The installer creates:
#   <VS_ROOT>/config/VS.config
#   <VS_ROOT>/config/VS.cfg
#   <VS_ROOT>/config/run.env
#   <VS_ROOT>/config/reads.txt
#   <VS_ROOT>/scripts/run_virusseeker.sh
#   <VS_ROOT>/scripts/validate_install.sh
#
# Database layout:
#   <VS_ROOT>/databases
#
# IMPORTANT:
#   vhunter_acc.db is built automatically from NCBI accession2taxid files.
#   taxa.sqlite is built automatically with ETE3 from the same NCBI taxdump
#   downloaded by this installer.
#
# Examples:
#
# Docker/local-sbatch:
#   ./install.sh \
#     --mode docker \
#     --root /export/apps/virusseeker \
#     --image virusseeker2:test \
#     --threads 30 \
#     --memory 64G \
#     --yes
#
# Native SLURM:
#   ./install.sh \
#     --mode slurm \
#     --root /export/apps/virusseeker \
#     --native-vs-home /path/to/VirusSeeker2.0 \
#     --conda-prefix /path/to/miniforge3 \
#     --slurm-partition genomics \
#     --threads 30 \
#     --memory 64G \
#     --yes

VS_EXECUTION_MODE="${VS_EXECUTION_MODE:-docker}"
VS_DB_PROFILE="${VS_DB_PROFILE:-vs}"
VS_ROOT="${VS_ROOT:-}"
VS_IMAGE="${VS_IMAGE:-virusseeker2:latest}"
THREADS="${THREADS:-$(nproc)}"
MEMORY="${MEMORY:-64G}"
MMSEQS_SPLIT_MEMORY="${MMSEQS_SPLIT_MEMORY:-100G}"
DONT_BUILD_NR_DIAMOND="${DONT_BUILD_NR_DIAMOND:-0}"

SLURM_PARTITION="${SLURM_PARTITION:-vs}"
NATIVE_VS_HOME="${NATIVE_VS_HOME:-}"
NATIVE_CONDA_PREFIX="${NATIVE_CONDA_PREFIX:-}"

YES=0

INSTALLER_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Dfam release defaults. Override these variables when publishing a newer release.
DFAM_BASE_URL="${DFAM_BASE_URL:-https://www.dfam.org/releases/current/families/FamDB}"
DFAM_ROOT_FILE="${DFAM_ROOT_FILE:-dfam40.0.h5.gz}"
DFAM_CURATED_FILE="${DFAM_CURATED_FILE:-dfam40.curated.consensus.0.h5.gz}"

usage() {
    cat <<'EOF'
VirusSeeker 2.0 installer

Required:
  --mode MODE              docker or slurm
  --db-profile PROFILE     test or vs
  --root PATH              Installation/runtime root

Docker mode:
  --image IMAGE            VirusSeeker Docker image

Native SLURM mode:
  --native-vs-home PATH    VirusSeeker2.0 Git checkout
  --conda-prefix PATH      Miniforge/Conda prefix used by native VirusSeeker
  --slurm-partition NAME   Site SLURM partition

Auxiliary databases:

General:
  --threads N              Threads/SLURM CPUs per task
  --memory SIZE            Memory, e.g. 64G
  --mmseqs-split-memory S  MMseqs split-memory limit, e.g. 100G
  --skip-nr-diamond        Do not build the full NR DIAMOND database
  --yes                    Confirm the very large database installation
  -h, --help               Show this help

Environment-variable equivalents:
  VS_EXECUTION_MODE
  VS_DB_PROFILE
  VS_ROOT
  VS_IMAGE
  THREADS
  MEMORY
  MMSEQS_SPLIT_MEMORY
  SLURM_PARTITION
  NATIVE_VS_HOME
  NATIVE_CONDA_PREFIX

Examples:

  Docker:
    ./install.sh \
      --mode docker \
      --root /export/apps/virusseeker \
      --image virusseeker2:test \
#     --threads 30 \
      --memory 64G \
      --yes

  Native SLURM:
    ./install.sh \
      --mode slurm \
      --root /export/apps/virusseeker \
      --native-vs-home "$PWD" \
      --conda-prefix /data/miniforge \
      --slurm-partition genomics \
#     --threads 30 \
      --memory 64G \
      --yes
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            [[ $# -ge 2 ]] || { echo "ERROR: --mode requires docker or slurm" >&2; exit 2; }
            VS_EXECUTION_MODE="$2"
            shift 2
            ;;
        --db-profile)
            [[ $# -ge 2 ]] || { echo "ERROR: --db-profile requires test or vs" >&2; exit 2; }
            VS_DB_PROFILE="$2"
            shift 2
            ;;
        --root)
            [[ $# -ge 2 ]] || { echo "ERROR: --root requires a path" >&2; exit 2; }
            VS_ROOT="$2"
            shift 2
            ;;
        --image)
            [[ $# -ge 2 ]] || { echo "ERROR: --image requires an image name" >&2; exit 2; }
            VS_IMAGE="$2"
            shift 2
            ;;
        --threads)
            [[ $# -ge 2 ]] || { echo "ERROR: --threads requires an integer" >&2; exit 2; }
            THREADS="$2"
            shift 2
            ;;
        --memory)
            [[ $# -ge 2 ]] || { echo "ERROR: --memory requires a value such as 64G" >&2; exit 2; }
            MEMORY="$2"
            shift 2
            ;;
        --mmseqs-split-memory)
            [[ $# -ge 2 ]] || { echo "ERROR: --mmseqs-split-memory requires a value" >&2; exit 2; }
            MMSEQS_SPLIT_MEMORY="$2"
            shift 2
            ;;
        --slurm-partition)
            [[ $# -ge 2 ]] || { echo "ERROR: --slurm-partition requires a value" >&2; exit 2; }
            SLURM_PARTITION="$2"
            shift 2
            ;;
        --native-vs-home)
            [[ $# -ge 2 ]] || { echo "ERROR: --native-vs-home requires a path" >&2; exit 2; }
            NATIVE_VS_HOME="$2"
            shift 2
            ;;
        --conda-prefix)
            [[ $# -ge 2 ]] || { echo "ERROR: --conda-prefix requires a path" >&2; exit 2; }
            NATIVE_CONDA_PREFIX="$2"
            shift 2
            ;;
        --skip-nr-diamond)
            DONT_BUILD_NR_DIAMOND=1
            shift
            ;;
        -y|--yes)
            YES=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$VS_EXECUTION_MODE" in
    docker|slurm) ;;
    *)
        echo "ERROR: --mode must be docker or slurm." >&2
        exit 2
        ;;
esac

case "$VS_DB_PROFILE" in
    test|vs) ;;
    *)
        echo "ERROR: --db-profile must be test or vs." >&2
        exit 2
        ;;
esac

[[ -n "$VS_ROOT" ]] || {
    echo "ERROR: choose an installation root with --root PATH." >&2
    exit 2
}

if [[ "$VS_ROOT" != /* ]]; then
    VS_ROOT="$(cd "$(dirname "$VS_ROOT")" 2>/dev/null && pwd)/$(basename "$VS_ROOT")"
fi
VS_ROOT="${VS_ROOT%/}"

# Avoid unsafe path rewriting in native mode.
[[ "$VS_ROOT" != *" "* ]] || {
    echo "ERROR: VS_ROOT cannot contain spaces." >&2
    exit 2
}

if [[ "$VS_EXECUTION_MODE" == "slurm" ]]; then
    [[ -n "$NATIVE_VS_HOME" ]] || {
        echo "ERROR: --native-vs-home is required in slurm mode." >&2
        exit 2
    }
    [[ -n "$NATIVE_CONDA_PREFIX" ]] || {
        echo "ERROR: --conda-prefix is required in slurm mode." >&2
        exit 2
    }
    NATIVE_VS_HOME="$(cd "$NATIVE_VS_HOME" && pwd)"
    NATIVE_CONDA_PREFIX="$(cd "$NATIVE_CONDA_PREFIX" && pwd)"
fi

[[ "$YES" -eq 1 ]] || {
    if [[ "$VS_DB_PROFILE" == "vs" ]]; then
        cat >&2 <<EOF
The VS database profile performs a fresh download/build of very large
biological databases and can require well over 1 TB including temporary data.

Run deliberately with --yes.
EOF
    else
        cat >&2 <<EOF
The test database profile builds compact format-compatible databases and
a compact public SARS-CoV-2 test data set.

Run deliberately with --yes.
EOF
    fi
    exit 2
}

parent_dir="$(dirname "$VS_ROOT")"
if [[ -e "$VS_ROOT" ]]; then
    [[ -d "$VS_ROOT" ]] || {
        echo "ERROR: VS_ROOT exists but is not a directory: $VS_ROOT" >&2
        exit 1
    }
    [[ -w "$VS_ROOT" ]] || {
        echo "ERROR: VS_ROOT is not writable by $(id -un): $VS_ROOT" >&2
        exit 1
    }
else
    [[ -d "$parent_dir" && -w "$parent_dir" ]] || {
        echo "ERROR: Cannot create $VS_ROOT as $(id -un)." >&2
        exit 1
    }
fi

if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
    command -v docker >/dev/null 2>&1 || {
        echo "ERROR: docker is required in docker mode." >&2
        exit 1
    }
    docker image inspect "$VS_IMAGE" >/dev/null 2>&1 || {
        echo "ERROR: Docker image is not available locally: $VS_IMAGE" >&2
        exit 1
    }
else
    command -v sbatch >/dev/null 2>&1 || {
        echo "ERROR: sbatch is required in slurm mode." >&2
        exit 1
    }
    [[ -r "$NATIVE_VS_HOME/queue.sh" ]] || {
        echo "ERROR: missing $NATIVE_VS_HOME/queue.sh" >&2
        exit 1
    }
    [[ -r "$NATIVE_VS_HOME/scripts/runVS.pl" ]] || {
        echo "ERROR: missing $NATIVE_VS_HOME/scripts/runVS.pl" >&2
        exit 1
    }
fi

INSTALL_UID="$(id -u)"
INSTALL_GID="$(id -g)"

DB_ROOT="$VS_ROOT/databases"
CONFIG_DIR="$VS_ROOT/config"
SCRIPT_DIR="$VS_ROOT/scripts"
SCRATCH_DIR="$VS_ROOT/scratch"
STATE_DIR="$VS_ROOT/state"
MARKER_DIR="$STATE_DIR/install-markers"
LOG_DIR="$STATE_DIR/install-logs"
RELEASE_DIR="$DB_ROOT/releases"
DATE_TAG="$(date +%Y%m%d)"

mkdir -p \
    "$VS_ROOT/input" \
    "$VS_ROOT/output" \
    "$CONFIG_DIR" \
    "$SCRIPT_DIR" \
    "$SCRATCH_DIR/mmseqs" \
    "$SCRATCH_DIR/diamond" \
    "$SCRATCH_DIR/build" \
    "$STATE_DIR" \
    "$MARKER_DIR" \
    "$LOG_DIR" \
    "$DB_ROOT/famdb" \
    "$DB_ROOT/ref" \
    "$RELEASE_DIR"

chown -R "$INSTALL_UID:$INSTALL_GID" "$VS_ROOT"
chmod -R u+rwX,go+rX "$VS_ROOT"

exec > >(tee -a "$LOG_DIR/install_${DATE_TAG}.log") 2>&1

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    log "ERROR: $*"
    exit 1
}

stage() {
    local name="$1"
    shift
    local marker="$MARKER_DIR/${name}.done"

    if [[ -f "$marker" ]]; then
        log "SKIP $name (marker exists)"
        return 0
    fi

    log "START $name"
    "$@"
    touch "$marker"
    chown "$INSTALL_UID:$INSTALL_GID" "$marker"
    log "DONE $name"
}

repair_test_profile_markers() {
    [[ "$VS_DB_PROFILE" == "test" ]] || return 0

    local test_db_marker="$MARKER_DIR/test_databases.done"
    local tax_marker="$MARKER_DIR/test_taxonomy.done"
    local vhunter_marker="$MARKER_DIR/vhunter_acc_db.done"

    # An older successful marker is not sufficient if required databases are
    # missing.  Force the deterministic tiny test DB stage to rebuild.
    if [[ -f "$test_db_marker" ]]; then
        if [[ ! -s "$DB_ROOT/core_nt_mmseqs/core_nt.dbtype" ]] || \
           [[ ! -s "$DB_ROOT/VirusDBNT/virus_nt.dbtype" ]] || \
           [[ ! -s "$DB_ROOT/VirusDBNR/virus_nr.dmnd" ]] || \
           [[ ! -s "$DB_ROOT/nr_dmnd/nr.dmnd" ]] || \
           [[ ! -s "$VS_ROOT/input/SRR22470068_R1.fastq.gz" ]] || \
           [[ ! -s "$VS_ROOT/input/SRR22470068_R2.fastq.gz" ]]; then
            log "RESET test_databases marker: required test database or VirusSeeker-compatible paired reads are missing"
            rm -f "$test_db_marker"
        fi
    fi

    if [[ -f "$CONFIG_DIR/reads.txt" ]] && \
       grep -Eq 'SRR22470068_[12]\.fastq\.gz' "$CONFIG_DIR/reads.txt"; then
        log "RESET reads_template marker: converting stale _1/_2 paired-read naming to _R1/_R2"
        rm -f "$MARKER_DIR/reads_template.done" "$CONFIG_DIR/reads.txt"
    fi


    if [[ -f "$tax_marker" ]]; then
        local ntmap="$DB_ROOT/taxdump/nucl_gb.accession2taxid"
        local protmap="$DB_ROOT/taxdump/prot.accession2taxid"
        if [[ ! -s "$ntmap" ]] || \
           [[ ! -s "$protmap" ]] || \
           ! grep -q $'NC_045512\tNC_045512.2\t2697049\t0' "$ntmap" || \
           ! grep -q $'NC_012920\tNC_012920.1\t9606\t0' "$ntmap" || \
           ! grep -q $'YP_009724389\tYP_009724389.1\t2697049\t0' "$protmap"; then
            log "RESET test_taxonomy/vhunter markers: public SARS-CoV-2 test mappings are missing or stale"
            rm -f "$tax_marker" "$vhunter_marker"
        fi
    fi

    # vhunter is tiny in test mode; if taxonomy was refreshed, rebuild it.
    if [[ ! -f "$tax_marker" ]]; then
        rm -f "$vhunter_marker"
    fi
}

vs_exec() {
    local command_text="$1"

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        docker run --rm -i \
            --platform linux/amd64 \
            --user "$INSTALL_UID:$INSTALL_GID" \
            --entrypoint /bin/bash \
            -e HOME=/tmp \
            -e TMPDIR=/scratch/build \
            -e FAMDB_DATA_DIR=/database/famdb \
            -v "$DB_ROOT:/database" \
            -v "$SCRATCH_DIR:/scratch" \
            -v "$VS_ROOT/input:/input" \
            -v "$INSTALLER_REPO_ROOT/installer:/installer:ro" \
            "$VS_IMAGE" \
            -lc "$command_text"
    else
        command_text="${command_text//\/database/$DB_ROOT}"
        command_text="${command_text//\/scratch/$SCRATCH_DIR}"
        command_text="${command_text//\/input/$VS_ROOT/input}"
        command_text="${command_text//\/installer/$INSTALLER_REPO_ROOT/installer}"
        HOME="${HOME:-/tmp}" \
        TMPDIR="$SCRATCH_DIR/build" \
        FAMDB_DATA_DIR="$DB_ROOT/famdb" \
        bash -lc "$command_text"
    fi
}

validate_toolchain() {
    vs_exec '
        set -Eeuo pipefail
        for command in \
            update_blastdb.pl blastdbcmd makeblastdb mmseqs diamond taxonkit \
            bowtie2-build pigz wget curl md5sum; do
            command -v "$command" >/dev/null || {
                echo "Missing command in image: $command" >&2
                exit 1
            }
        done
        [[ -x /opt/conda/envs/vs/bin/fasterq-dump ]] || {
            echo "Missing /opt/conda/envs/vs/bin/fasterq-dump in image" >&2
            exit 1
        }
        /opt/conda/envs/vs/bin/python -c "from Bio import SeqIO"
        /opt/conda/envs/vs/bin/python -c "from ete3 import NCBITaxa"
        echo "Container database toolchain, Biopython, and ETE3 are available."
    '
}


download_test_taxonomy() {
    local release_host="$RELEASE_DIR/taxdump_test_${DATE_TAG}"
    local release_container="/database/releases/taxdump_test_${DATE_TAG}"

    mkdir -p "$release_host"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$release_host"

    vs_exec "/installer/build_test_taxonomy.sh '$release_container'"

    ln -sfn "releases/taxdump_test_${DATE_TAG}" "$DB_ROOT/taxdump"
    ln -sfn "taxdump/virus.taxid.txt" "$DB_ROOT/virus.taxid.txt"
}

build_test_databases_and_reads() {
    vs_exec '
        /installer/download_real_test_assets.sh
        /installer/build_test_databases.sh
    '
}

build_virus_family_genomes() {
    local out="$DB_ROOT/virus_family_genomes.tsv"

    if [[ "$VS_DB_PROFILE" == "test" ]]; then
        local sars_fasta="$DB_ROOT/test_source/NC_045512.2.fa"
        [[ -s "$sars_fasta" ]] || die "Missing real SARS-CoV-2 reference: $sars_fasta"
        local genome_size
        genome_size="$(awk '!/^>/ {gsub(/[[:space:]]/,""); n+=length($0)} END {print n+0}' "$sars_fasta")"
        [[ "$genome_size" -gt 0 ]] || die "Could not determine SARS-CoV-2 genome size"
        printf 'fam\tsize\nCoronaviridae\t%s\n' "$genome_size" > "$out"
        chown "$INSTALL_UID:$INSTALL_GID" "$out"
        log "Created real-data test virus-family genome-size table: $out ($genome_size bp)"
        return 0
    fi

    if [[ -s "$out" ]]; then
        log "Using existing production virus-family genome-size table: $out"
        return 0
    fi

    die "Production profile requires $out (two-column TSV: fam<TAB>size, derived from ICTV Virus Properties)"
}

write_test_expected_results() {
    mkdir -p "$VS_ROOT/test_data"
    install -m 0644 \
        "$INSTALLER_REPO_ROOT/installer/templates/expected_results.txt" \
        "$VS_ROOT/test_data/expected_results.txt"
    cp "$VS_ROOT/input/README.test-data.txt" "$VS_ROOT/test_data/README.txt"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$VS_ROOT/test_data"
}

download_blast_db() {
    local database="$1"
    local release="/database/releases/${database}_${DATE_TAG}"

    mkdir -p "$DB_ROOT/releases/${database}_${DATE_TAG}"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$DB_ROOT/releases/${database}_${DATE_TAG}"

    vs_exec "
        set -Eeuo pipefail
        cd '$release'
        update_blastdb.pl --decompress '$database'
        test -e '$release/${database}.pal' \
          -o -e '$release/${database}.nal' \
          -o -e '$release/${database}.00.psq' \
          -o -e '$release/${database}.00.nsq' \
          -o -e '$release/${database}.psq' \
          -o -e '$release/${database}.nsq'
    "

    ln -sfn "releases/${database}_${DATE_TAG}" "$DB_ROOT/$database"
}

download_nr() {
    download_blast_db nr
}

download_nt() {
    download_blast_db nt
}

download_core_nt() {
    download_blast_db core_nt
}

download_famdb_file() {
    local filename="$1"
    local destination="$DB_ROOT/famdb/$filename"
    local checksum="$destination.md5"

    vs_exec "
        set -Eeuo pipefail
        cd /database/famdb
        wget -c -O '$filename' '$DFAM_BASE_URL/$filename'
        wget -c -O '$filename.md5' '$DFAM_BASE_URL/$filename.md5'
        expected=\$(awk '{print \$1}' '$filename.md5')
        actual=\$(md5sum '$filename' | awk '{print \$1}')
        [[ \"\$expected\" == \"\$actual\" ]] || {
            echo 'Checksum mismatch for $filename' >&2
            exit 1
        }
        pigz -dkf '$filename'
    "

    [[ -s "${destination%.gz}" ]] ||
        die "FamDB decompression failed: ${destination%.gz}"
}

download_famdb() {
    download_famdb_file "$DFAM_ROOT_FILE"
    download_famdb_file "$DFAM_CURATED_FILE"
}

download_human_reference() {
    local release_host="$RELEASE_DIR/ref_GRCh38.p14_${DATE_TAG}"
    local release_container="/database/releases/ref_GRCh38.p14_${DATE_TAG}"
    local filename="GCF_000001405.40_GRCh38.p14_genomic.fna"
    local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"

    mkdir -p "$release_host"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$release_host"

    vs_exec "
        set -Eeuo pipefail
        cd '$release_container'
        wget -c -O '$filename.gz' '$url'
        pigz -dkf '$filename.gz'
        test -s '$filename'
        bowtie2-build --threads '$THREADS' '$filename' GRCh38.p14
    "

    rm -rf "$DB_ROOT/ref"
    ln -sfn "releases/ref_GRCh38.p14_${DATE_TAG}" "$DB_ROOT/ref"
}

download_taxonomy() {
    local release_host="$RELEASE_DIR/taxdump_${DATE_TAG}"
    local release_container="/database/releases/taxdump_${DATE_TAG}"

    mkdir -p "$release_host"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$release_host"

    vs_exec "
        set -Eeuo pipefail
        cd '$release_container'
        rm -f taxdump.tar.gz.part

        if [[ -s taxdump.tar.gz ]] && gzip -t taxdump.tar.gz 2>/dev/null; then
            echo 'Existing taxdump.tar.gz passed gzip validation; reusing it.'
        else
            rm -f taxdump.tar.gz
            wget --tries=5 --timeout=30 \
              -O taxdump.tar.gz.part \
              https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

            gzip -t taxdump.tar.gz.part
            mv -f taxdump.tar.gz.part taxdump.tar.gz
        fi

        gzip -t taxdump.tar.gz
        tar -xzf taxdump.tar.gz

        wget -c -O nucl_gb.accession2taxid.gz \
          https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        wget -c -O prot.accession2taxid.gz \
          https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

        pigz -dkf nucl_gb.accession2taxid.gz
        pigz -dkf prot.accession2taxid.gz

        taxonkit list \
          --data-dir '$release_container' \
          --ids 10239 \
          --indent '' \
          > virus.taxid.txt

        test -s virus.taxid.txt
    "

    ln -sfn "releases/taxdump_${DATE_TAG}" "$DB_ROOT/taxdump"
    ln -sfn "taxdump/virus.taxid.txt" "$DB_ROOT/virus.taxid.txt"
}


build_vhunter_acc_db() {
    local nucl_file="$DB_ROOT/taxdump/nucl_gb.accession2taxid"
    local prot_file="$DB_ROOT/taxdump/prot.accession2taxid"
    local output="$DB_ROOT/vhunter_acc.db"
    local tmp_output="$DB_ROOT/vhunter_acc.db.tmp"

    [[ -s "$nucl_file" ]] || die "Missing $nucl_file"
    [[ -s "$prot_file" ]] || die "Missing $prot_file"
    rm -f "$tmp_output"

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        docker run --rm \
            --platform linux/amd64 \
            --user "$INSTALL_UID:$INSTALL_GID" \
            --entrypoint /opt/conda/envs/vs/bin/python \
            -v "$DB_ROOT:/database" \
            -v "$INSTALLER_REPO_ROOT/installer:/installer:ro" \
            "$VS_IMAGE" \
            /installer/build_vhunter_acc.py \
                --nucl /database/taxdump/nucl_gb.accession2taxid \
                --prot /database/taxdump/prot.accession2taxid \
                --output /database/vhunter_acc.db.tmp
    else
        local python_bin="$NATIVE_CONDA_PREFIX/bin/python"
        [[ -x "$python_bin" ]] || python_bin="$(command -v python3 || command -v python)"
        [[ -n "$python_bin" ]] || die "Python is required to build vhunter_acc.db"

        "$python_bin" "$INSTALLER_REPO_ROOT/installer/build_vhunter_acc.py" \
            --nucl "$nucl_file" \
            --prot "$prot_file" \
            --output "$tmp_output"
    fi

    [[ -s "$tmp_output" ]] || die "Failed to build $tmp_output"
    mv -f "$tmp_output" "$output"
    chmod 0644 "$output"

    if [[ "$VS_DB_PROFILE" == "test" ]]; then
        if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
            docker run --rm \
                --platform linux/amd64 \
                --entrypoint /opt/conda/envs/vs/bin/python \
                -v "$DB_ROOT:/database:ro" \
                -v "$INSTALLER_REPO_ROOT/installer:/installer:ro" \
                "$VS_IMAGE" \
                /installer/validate_sars_mappings.py /database/vhunter_acc.db
        else
            local python_bin="$NATIVE_CONDA_PREFIX/bin/python"
            [[ -x "$python_bin" ]] || python_bin="$(command -v python3 || command -v python)"
            "$python_bin" "$INSTALLER_REPO_ROOT/installer/validate_sars_mappings.py" "$output"
        fi
    fi

    log "Built $output"
}

build_taxa_sqlite() {
    local taxdump="$DB_ROOT/taxdump/taxdump.tar.gz"
    local output="$DB_ROOT/taxa.sqlite"
    local traverse="$DB_ROOT/taxa.sqlite.traverse.pkl"

    [[ -s "$taxdump" ]] || die "Missing NCBI taxonomy archive: $taxdump"
    rm -f "$output" "$traverse"

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        docker run --rm \
            --platform linux/amd64 \
            --user "$INSTALL_UID:$INSTALL_GID" \
            --entrypoint /opt/conda/envs/vs/bin/python \
            -e PYTHONNOUSERSITE=1 \
            -v "$DB_ROOT:/database" \
            -v "$INSTALLER_REPO_ROOT/installer:/installer:ro" \
            "$VS_IMAGE" \
            /installer/build_taxa_sqlite.py \
                --taxdump /database/taxdump/taxdump.tar.gz \
                --output /database/taxa.sqlite
    else
        local python_bin="$NATIVE_CONDA_PREFIX/bin/python"
        [[ -x "$python_bin" ]] || python_bin="$(command -v python3 || true)"
        [[ -n "$python_bin" && -x "$python_bin" ]] || die "Python is required to build taxa.sqlite"

        PYTHONNOUSERSITE=1 "$python_bin" \
            "$INSTALLER_REPO_ROOT/installer/build_taxa_sqlite.py" \
            --taxdump "$taxdump" \
            --output "$output"
    fi

    [[ -s "$output" ]] || die "Failed to create $output"
    chmod 0644 "$output"
    [[ ! -e "$traverse" ]] || chmod 0644 "$traverse"
    log "Built ETE3 taxonomy database: $output"
}

build_virus_nr() {
    local release_host="$RELEASE_DIR/VirusDBNR_${DATE_TAG}"
    local release_container="/database/releases/VirusDBNR_${DATE_TAG}"

    mkdir -p "$release_host" "$SCRATCH_DIR/mmseqs/virus_nr"
    chown -R "$INSTALL_UID:$INSTALL_GID" \
        "$release_host" "$SCRATCH_DIR/mmseqs/virus_nr"

    vs_exec "
        set -Eeuo pipefail

        blastdbcmd \
          -taxidlist /database/virus.taxid.txt \
          -dbtype prot \
          -db /database/nr/nr \
          > '$release_container/virus_nr.fasta'

        test -s '$release_container/virus_nr.fasta'

        mmseqs easy-linclust \
          '$release_container/virus_nr.fasta' \
          '$release_container/virus_nr.clustr98_98' \
          /scratch/mmseqs/virus_nr \
          --min-seq-id 0.98 \
          -c 0.98 \
          --threads '$THREADS'

        test -s '$release_container/virus_nr.clustr98_98_rep_seq.fasta'

        diamond makedb \
          --in '$release_container/virus_nr.clustr98_98_rep_seq.fasta' \
          --db '$release_container/virus_nr'

        test -s '$release_container/virus_nr.dmnd'
    "

    ln -sfn "releases/VirusDBNR_${DATE_TAG}" "$DB_ROOT/VirusDBNR"
}

build_virus_nt() {
    local release_host="$RELEASE_DIR/VirusDBNT_${DATE_TAG}"
    local release_container="/database/releases/VirusDBNT_${DATE_TAG}"

    mkdir -p "$release_host" "$SCRATCH_DIR/mmseqs/virus_nt"
    chown -R "$INSTALL_UID:$INSTALL_GID" \
        "$release_host" "$SCRATCH_DIR/mmseqs/virus_nt"

    vs_exec "
        set -Eeuo pipefail

        blastdbcmd \
          -taxidlist /database/virus.taxid.txt \
          -dbtype nucl \
          -db /database/core_nt/core_nt \
          > '$release_container/virus_nt.fasta'

        test -s '$release_container/virus_nt.fasta'

        mmseqs easy-linclust \
          '$release_container/virus_nt.fasta' \
          '$release_container/virus_nt.clustr98_98' \
          /scratch/mmseqs/virus_nt \
          --split-memory-limit '$MMSEQS_SPLIT_MEMORY' \
          --min-seq-id 0.98 \
          -c 0.98 \
          --threads '$THREADS'

        test -s '$release_container/virus_nt.clustr98_98_rep_seq.fasta'

        mmseqs createdb \
          '$release_container/virus_nt.clustr98_98_rep_seq.fasta' \
          '$release_container/virus_nt'

        test -s '$release_container/virus_nt.dbtype'
    "

    ln -sfn "releases/VirusDBNT_${DATE_TAG}" "$DB_ROOT/VirusDBNT"

    # runVS.pl expects ref/ref_viruses_rep_genomes beneath the configured ref path.
    cp -fL       "$release_host/virus_nt.clustr98_98_rep_seq.fasta"       "$DB_ROOT/ref/ref_viruses_rep_genomes"

    vs_exec "
        set -Eeuo pipefail
        bowtie2-build           /database/ref/ref_viruses_rep_genomes           /database/ref/ref_viruses_rep_genomes
    "
}

build_core_nt_mmseqs() {
    local release_host="$RELEASE_DIR/core_nt_mmseqs_${DATE_TAG}"
    local release_container="/database/releases/core_nt_mmseqs_${DATE_TAG}"

    mkdir -p "$release_host"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$release_host"

    vs_exec "
        set -Eeuo pipefail

        blastdbcmd \
          -dbtype nucl \
          -db /database/core_nt/core_nt \
          -entry all \
          | pigz -p '$THREADS' \
          > '$release_container/core_nt.fasta.gz'

        test -s '$release_container/core_nt.fasta.gz'

        mmseqs createdb \
          '$release_container/core_nt.fasta.gz' \
          '$release_container/core_nt'

        test -s '$release_container/core_nt.dbtype'
    "

    ln -sfn "releases/core_nt_mmseqs_${DATE_TAG}" "$DB_ROOT/core_nt_mmseqs"
}

build_nr_diamond() {
    local release_host="$RELEASE_DIR/nr_dmnd_${DATE_TAG}"
    local release_container="/database/releases/nr_dmnd_${DATE_TAG}"

    mkdir -p "$release_host"
    chown -R "$INSTALL_UID:$INSTALL_GID" "$release_host"

    vs_exec "
        set -Eeuo pipefail

        blastdbcmd \
          -dbtype prot \
          -db /database/nr/nr \
          -entry all \
          | pigz -p '$THREADS' \
          > '$release_container/nr.fasta.gz'

        test -s '$release_container/nr.fasta.gz'

        diamond makedb \
          --in '$release_container/nr.fasta.gz' \
          --db '$release_container/nr'

        test -s '$release_container/nr.dmnd'
    "

    ln -sfn "releases/nr_dmnd_${DATE_TAG}" "$DB_ROOT/nr_dmnd"
}

write_vs_config() {
    local runtime_db runtime_scratch runtime_conda runtime_bin partition

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        runtime_db="/database"
        runtime_scratch="/scratch"
        runtime_conda="/opt/conda"
        runtime_bin="/opt/virusseeker/third_party/"
        partition="local"
    else
        runtime_db="$DB_ROOT"
        runtime_scratch="$SCRATCH_DIR"
        runtime_conda="$NATIVE_CONDA_PREFIX"
        runtime_bin="$NATIVE_VS_HOME/third_party/"
        partition="$SLURM_PARTITION"
    fi

    # Main VirusSeeker configuration used by runVS.pl.
    cat > "$CONFIG_DIR/VS.config" <<EOF
slurm_cpus_per_task=$THREADS
slurm_mem="$MEMORY"
slurm_partition="$partition"

diamond_args="--block-size 10 --iterate faster --index-chunks 1 --tmpdir $runtime_scratch/diamond"

bin_path="$runtime_bin"
conda_prefix="$runtime_conda"

ref="$runtime_db/ref"
nr="$runtime_db/nr_dmnd/nr.dmnd"
mmseqsnt="$runtime_db/core_nt_mmseqs/core_nt"
mmseqs_tmpdir="$runtime_scratch/mmseqs"
virus_nr="$runtime_db/VirusDBNR/virus_nr.dmnd"
virus_nt="$runtime_db/VirusDBNT/virus_nt"

# Python helper scripts in this VirusSeeker source read these FLAT keys
# directly from VS.config (not from VS.cfg).
# Do not quote these values: the parser does not strip surrounding quotes.
vhunter=$runtime_db/vhunter_acc.db
ncbi_taxdb=$runtime_db/taxa.sqlite
taxdump=$runtime_db/taxdump
virus_family_genomes=$runtime_db/virus_family_genomes.tsv
EOF

    cat > "$CONFIG_DIR/VS.cfg" <<EOF
[slurm]
slurm_cpus_per_task = $THREADS
slurm_mem = $MEMORY
slurm_partition = $partition

[diamond]
diamond_args = --block-size 10 --iterate faster --index-chunks 1 --tmpdir $runtime_scratch/diamond

[paths]
conda_prefix = $runtime_conda
db_dir = $runtime_db
taxdump = %(db_dir)s/taxdump
vhunter = %(db_dir)s/vhunter_acc.db
ncbi_taxadb = %(db_dir)s/taxa.sqlite
EOF

    chown "$INSTALL_UID:$INSTALL_GID" \
        "$CONFIG_DIR/VS.config" "$CONFIG_DIR/VS.cfg"
    chmod 0644 "$CONFIG_DIR/VS.config" "$CONFIG_DIR/VS.cfg"
}


validate_runtime_source() {
    local source_runvs="$INSTALLER_REPO_ROOT/scripts/runVS.pl"

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        docker run --rm \
            --platform linux/amd64 \
            --entrypoint /bin/bash \
            "$VS_IMAGE" \
            -lc '
                set -Eeuo pipefail
                f=/opt/virusseeker/scripts/runVS.pl
                test -s "$f"
                grep -q "bbsplit.sh path=.*bbmap_viral_ref_index" "$f"
                grep -q "bbmap.sh path=.*bbmap_host_index" "$f"
                grep -E -q "if \[ -s .*IN1.*\] && \[ -s .*IN2.*\]" "$f"
                grep -q "No paired reads remain after assembly; skipping PEAR cleanly." "$f"
                /opt/conda/envs/vs/bin/perl -c "$f" >/dev/null

                p=/opt/virusseeker/scripts/parse_blast_Virus_noBP.py
                test -s "$p"
                grep -q "No BLAST/MMseqs hits found; treating all input sequences as unassigned." "$p"
                /opt/conda/envs/vs/bin/python -m py_compile "$p"

                echo "Baked-in VirusSeeker source fixes: OK"
            '
    else
        source_runvs="$NATIVE_VS_HOME/scripts/runVS.pl"
        [[ -s "$source_runvs" ]] || die "Missing VirusSeeker source: $source_runvs"

        grep -q 'bbsplit.sh path=.*bbmap_viral_ref_index' "$source_runvs" || \
            die "Native runVS.pl is missing the BBSplitter writable-index fix"
        grep -q 'bbmap.sh path=.*bbmap_host_index' "$source_runvs" || \
            die "Native runVS.pl is missing the host BBMap writable-index fix"
        grep -E -q 'if \[ -s .*IN1.*\] && \[ -s .*IN2.*\]' "$source_runvs" || \
            die "Native runVS.pl is missing the metaSPAdes input guard"
        grep -q 'No paired reads remain after assembly; skipping PEAR cleanly.' "$source_runvs" || \
            die "Native runVS.pl is missing the empty-PEAR guard"

        local parser="$NATIVE_VS_HOME/scripts/parse_blast_Virus_noBP.py"
        [[ -s "$parser" ]] || die "Missing native parser: $parser"
        grep -q 'No BLAST/MMseqs hits found; treating all input sequences as unassigned.' "$parser" || \
            die "Native parser is missing empty-search handling"
    fi

    log "Validated application fixes from source/image; installer performed no source patching"
}

write_run_env() {
    local runtime_vhunter runtime_taxadb runtime_taxdump runtime_conda partition
    local run_name read_type host_db repeat_library

    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        runtime_vhunter="/database/vhunter_acc.db"
        runtime_taxadb="/database/taxa.sqlite"
        runtime_taxdump="/database/taxdump"
        runtime_conda="/opt/conda"
        partition="local"
    else
        runtime_vhunter="$DB_ROOT/vhunter_acc.db"
        runtime_taxadb="$DB_ROOT/taxa.sqlite"
        runtime_taxdump="$DB_ROOT/taxdump"
        runtime_conda="$NATIVE_CONDA_PREFIX"
        partition="$SLURM_PARTITION"
    fi

    if [[ "$VS_DB_PROFILE" == "test" ]]; then
        run_name="toy_test"
        read_type="s"

        # Full toy/test workflow:
        #   - non-X host reference => stage 3 host removal runs
        #   - Discovery mode (no -v) => stage 4 assembly runs
        if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
            host_db="/database/ref/GRCh38.p14"
            repeat_library="/database/test/repeat_library.fa"
        else
            host_db="$DB_ROOT/ref/GRCh38.p14"
            repeat_library="$DB_ROOT/test/repeat_library.fa"
        fi
    else
        run_name="real_run_001"
        read_type="s"
        host_db="X"
        repeat_library=""
    fi

    cat > "$CONFIG_DIR/run.env" <<EOF
# Generated by VirusSeeker install.sh
export VS_EXECUTION_MODE="$VS_EXECUTION_MODE"
export VS_DB_PROFILE="$VS_DB_PROFILE"

# Common installation paths
export VS_ROOT="$VS_ROOT"
export VS_INPUT_DIR="\$VS_ROOT/input"
export VS_OUTPUT_DIR="\$VS_ROOT/output"
export VS_CONFIG_DIR="\$VS_ROOT/config"
export VS_SCRATCH_DIR="\$VS_ROOT/scratch"
export VS_STATE_DIR="\$VS_ROOT/state"
export VS_DATABASE_DIR="\$VS_ROOT/databases"

export VS_READS_FILE="\$VS_CONFIG_DIR/reads.txt"
export VS_CONFIG_FILE="\$VS_CONFIG_DIR/VS.config"
export VS_CFG_FILE="\$VS_CONFIG_DIR/VS.cfg"

# Some legacy helper scripts resolve VS.config relative to their own script
# directory rather than from the parent scripts/ directory.
export VS_READ_COUNTER_CONFIG_RUNTIME="/opt/virusseeker/scripts/VS_Read_Counter/VS.config"

# Fixed runtime locations expected by VirusSeeker source code.
# parse_blast_Virus_noBP.py resolves VS.cfg next to the Python script.
if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
    export VS_CONFIG_RUNTIME="/opt/virusseeker/scripts/VS.config"
    export VS_CFG_RUNTIME="/opt/virusseeker/scripts/VS.cfg"
else
    export VS_CONFIG_RUNTIME="$NATIVE_VS_HOME/scripts/VS.config"
    export VS_CFG_RUNTIME="$NATIVE_VS_HOME/scripts/VS.cfg"
fi

# Runtime/database settings
export VS_VHUNTER_HOST="\$VS_DATABASE_DIR/vhunter_acc.db"
export VS_NCBI_TAXADB_HOST="\$VS_DATABASE_DIR/taxa.sqlite"
export VS_TAXDUMP_HOST="\$VS_DATABASE_DIR/taxdump"

export VS_VHUNTER_RUNTIME="$runtime_vhunter"
export VS_NCBI_TAXADB_RUNTIME="$runtime_taxadb"
export VS_TAXDUMP_RUNTIME="$runtime_taxdump"

export VS_CONDA_PREFIX="$runtime_conda"
export VS_SLURM_PARTITION="$partition"
export VS_REPEATMASKER_LIBRARY="$repeat_library"

# Run settings -- edit these for each run.
export VS_RUN_NAME="$run_name"
export VS_HOST_DB="$host_db"
export VS_READ_TYPE="$read_type"
export VS_ANALYSIS_MODE="m"
export VS_MEMORY="$MEMORY"
export VS_CPUS="$THREADS"

# Test profile deliberately uses Discovery mode so assembly executes.
# Production keeps the prior installer default (Virome mode) unless edited.
if [[ "$VS_DB_PROFILE" == "test" ]]; then
    export VS_EXTRA_FLAGS=""
else
    export VS_EXTRA_FLAGS="-v"
fi

# Docker mode
export VS_IMAGE="$VS_IMAGE"
export VS_CONTAINER="virusseeker-\${VS_RUN_NAME}"

# Native SLURM mode
export VS_NATIVE_VS_HOME="$NATIVE_VS_HOME"
EOF

    chown "$INSTALL_UID:$INSTALL_GID" "$CONFIG_DIR/run.env"
    chmod 0644 "$CONFIG_DIR/run.env"
}

write_reads_template() {
    if [[ -e "$CONFIG_DIR/reads.txt" ]]; then
        return 0
    fi

    if [[ "$VS_DB_PROFILE" == "test" ]]; then
        if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
            cat > "$CONFIG_DIR/reads.txt" <<'EOF'
/input/SRR22470068_R1.fastq.gz
/input/SRR22470068_R2.fastq.gz
EOF
        else
            cat > "$CONFIG_DIR/reads.txt" <<EOF
$VS_ROOT/input/SRR22470068_R1.fastq.gz
$VS_ROOT/input/SRR22470068_R2.fastq.gz
EOF
        fi
    elif [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        cat > "$CONFIG_DIR/reads.txt" <<'EOF'
/input/sample_R1.fastq.gz
/input/sample_R2.fastq.gz
EOF
    else
        cat > "$CONFIG_DIR/reads.txt" <<EOF
$VS_ROOT/input/sample_R1.fastq.gz
$VS_ROOT/input/sample_R2.fastq.gz
EOF
    fi

    chown "$INSTALL_UID:$INSTALL_GID" "$CONFIG_DIR/reads.txt"
    chmod 0644 "$CONFIG_DIR/reads.txt"
}

write_run_script() {
    install -m 0755 \
        "$INSTALLER_REPO_ROOT/installer/templates/run_virusseeker.sh" \
        "$SCRIPT_DIR/run_virusseeker.sh"
    chown "$INSTALL_UID:$INSTALL_GID" "$SCRIPT_DIR/run_virusseeker.sh"
}

write_validate_script() {
    install -m 0755 \
        "$INSTALLER_REPO_ROOT/installer/templates/validate_install.sh" \
        "$SCRIPT_DIR/validate_install.sh"
    chown "$INSTALL_UID:$INSTALL_GID" "$SCRIPT_DIR/validate_install.sh"
}

write_database_layout() {
    install -m 0644 \
        "$INSTALLER_REPO_ROOT/installer/templates/database_README.md" \
        "$DB_ROOT/README.md"
    chown "$INSTALL_UID:$INSTALL_GID" "$DB_ROOT/README.md"
}

main() {
    log "VS_ROOT=$VS_ROOT"
    log "VS_IMAGE=$VS_IMAGE"
    log "VS_DB_PROFILE=$VS_DB_PROFILE"
    log "THREADS=$THREADS"
    log "MEMORY=$MEMORY"

    stage validate_toolchain validate_toolchain
    repair_test_profile_markers

    if [[ "$VS_DB_PROFILE" == "test" ]]; then
        stage test_taxonomy download_test_taxonomy
        stage test_databases build_test_databases_and_reads
        stage vhunter_acc_db build_vhunter_acc_db
        stage taxa_sqlite build_taxa_sqlite
        stage virus_family_genomes build_virus_family_genomes
        stage test_expected_results write_test_expected_results
    else
        stage famdb download_famdb
        stage nr download_nr
        stage nt download_nt
        stage core_nt download_core_nt
        stage human_reference download_human_reference
        stage taxonomy download_taxonomy
        stage vhunter_acc_db build_vhunter_acc_db
        stage taxa_sqlite build_taxa_sqlite
        stage virus_family_genomes build_virus_family_genomes
        stage virus_nr build_virus_nr
        stage virus_nt build_virus_nt
        stage core_nt_mmseqs build_core_nt_mmseqs

        if [[ "$DONT_BUILD_NR_DIAMOND" == "1" ]]; then
            log "SKIP nr_diamond because DONT_BUILD_NR_DIAMOND=1"
        else
            stage nr_diamond build_nr_diamond
        fi
    fi

    stage vs_config write_vs_config
    stage runtime_source validate_runtime_source
    stage run_env write_run_env
    stage reads_template write_reads_template
    stage run_script write_run_script
    stage validate_script write_validate_script
    stage database_layout write_database_layout

    chown -R "$INSTALL_UID:$INSTALL_GID" "$VS_ROOT"
    chmod -R u+rwX,go+rX "$VS_ROOT"

    cat <<EOF

VirusSeeker installation completed.

Execution mode:
  $VS_EXECUTION_MODE

Database profile:
  $VS_DB_PROFILE

Host root:
  $VS_ROOT

Configuration:
  $VS_ROOT/config/VS.config
  $VS_ROOT/config/VS.cfg
  $VS_ROOT/config/run.env
  $VS_ROOT/config/reads.txt

Validate:
  $VS_ROOT/scripts/validate_install.sh

Run:
  $VS_ROOT/scripts/run_virusseeker.sh $VS_ROOT/config/run.env
EOF
}

main "$@"
