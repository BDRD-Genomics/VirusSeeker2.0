#!/usr/bin/env bash
set -Eeuo pipefail

ENV_FILE="${1:-}"
if [[ -z "$ENV_FILE" ]]; then
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    ENV_FILE="$(cd "$script_dir/.." && pwd)/config/run.env"
fi

[[ -r "$ENV_FILE" ]] || {
    echo "ERROR: cannot read $ENV_FILE" >&2
    exit 1
}

source "$ENV_FILE"

required=(
    VS_EXECUTION_MODE VS_DB_PROFILE VS_ROOT VS_INPUT_DIR VS_OUTPUT_DIR VS_CONFIG_DIR
    VS_SCRATCH_DIR VS_STATE_DIR VS_DATABASE_DIR VS_READS_FILE
    VS_CONFIG_FILE VS_CFG_FILE VS_CONFIG_RUNTIME VS_CFG_RUNTIME
    VS_READ_COUNTER_CONFIG_RUNTIME
    VS_RUN_NAME VS_HOST_DB VS_READ_TYPE
    VS_ANALYSIS_MODE VS_MEMORY VS_CPUS
)

for name in "${required[@]}"; do
    [[ -n "${!name:-}" ]] || {
        echo "ERROR: missing variable $name" >&2
        exit 1
    }
done

[[ -s "$VS_VHUNTER_HOST" ]] || {
    echo "ERROR: missing $VS_VHUNTER_HOST" >&2
    exit 1
}
[[ -s "$VS_NCBI_TAXADB_HOST" ]] || {
    echo "ERROR: missing $VS_NCBI_TAXADB_HOST" >&2
    exit 1
}

[[ -s "$VS_CONFIG_FILE" ]] || {
    echo "ERROR: missing $VS_CONFIG_FILE" >&2
    exit 1
}
[[ -s "$VS_CFG_FILE" ]] || {
    echo "ERROR: missing $VS_CFG_FILE" >&2
    exit 1
}
# Validate BOTH config formats before starting any work.
# The current parse_blast_Virus_noBP.py reads flat key=value pairs from
# scripts/VS.config and specifically requires "vhunter" and "ncbi_taxdb".
python3 - "$VS_CONFIG_FILE" "$VS_CFG_FILE" <<'PY_CFG_CHECK'
import configparser
import sys

vs_config_path, vs_cfg_path = sys.argv[1], sys.argv[2]

flat = {}
with open(vs_config_path, "r", encoding="utf-8") as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        flat[key.strip()] = value.strip()

required_flat = ("vhunter", "ncbi_taxdb", "virus_family_genomes")
missing_flat = [key for key in required_flat if not flat.get(key, "")]
if missing_flat:
    raise SystemExit(
        "ERROR: VS.config missing required flat key(s): "
        + ", ".join(missing_flat)
    )

# These two paths must not contain shell-style quote characters because
# parse_blast_Virus_noBP.py uses the values literally.
for key in required_flat:
    value = flat[key]
    if value[:1] in ("'", '"') or value[-1:] in ("'", '"'):
        raise SystemExit(
            f"ERROR: VS.config {key} must be unquoted for Python parser: {value}"
        )

cfg = configparser.ConfigParser()
if not cfg.read(vs_cfg_path):
    raise SystemExit(f"ERROR: could not read VS.cfg: {vs_cfg_path}")
if "paths" not in cfg:
    raise SystemExit(f"ERROR: VS.cfg is missing [paths]: {vs_cfg_path}")

required_ini = ("conda_prefix", "db_dir", "taxdump", "vhunter", "ncbi_taxadb")
missing_ini = [key for key in required_ini if not cfg["paths"].get(key, "").strip()]
if missing_ini:
    raise SystemExit(
        "ERROR: VS.cfg [paths] missing required key(s): "
        + ", ".join(missing_ini)
    )

print("VS.config flat-key preflight OK:", vs_config_path)
print("  vhunter              =", flat["vhunter"])
print("  ncbi_taxdb           =", flat["ncbi_taxdb"])
print("  virus_family_genomes =", flat["virus_family_genomes"])
print("VS.cfg INI preflight OK:", vs_cfg_path)
PY_CFG_CHECK

if [[ "$VS_DB_PROFILE" == "test" ]]; then
    if [[ "$VS_EXECUTION_MODE" == "docker" ]]; then
        [[ "$VS_HOST_DB" == "/database/ref/GRCh38.p14" ]] || {
            echo "ERROR: test profile must use /database/ref/GRCh38.p14 for host removal" >&2
            exit 1
        }
    else
        [[ -s "$VS_HOST_DB" ]] || {
            echo "ERROR: test-profile host reference is missing: $VS_HOST_DB" >&2
            exit 1
        }
    fi

    [[ -z "${VS_EXTRA_FLAGS:-}" ]] || {
        echo "ERROR: test profile must run Discovery mode; remove -v from VS_EXTRA_FLAGS" >&2
        exit 1
    }
fi

RUN_DIR="$VS_OUTPUT_DIR/$VS_RUN_NAME"
[[ ! -e "$RUN_DIR" ]] || {
    echo "ERROR: output already exists: $RUN_DIR" >&2
    exit 1
}

mkdir -p \
    "$VS_SCRATCH_DIR/mmseqs" \
    "$VS_SCRATCH_DIR/diamond" \
    "$VS_STATE_DIR"

case "$VS_EXECUTION_MODE" in
    docker)
        command -v docker >/dev/null 2>&1 || {
            echo "ERROR: docker is not available." >&2
            exit 1
        }

        while IFS= read -r container_path; do
            container_path="${container_path%$'\r'}"
            [[ -z "$container_path" ]] && continue
            case "$container_path" in
                /input/*) ;;
                *)
                    echo "ERROR: Docker-mode reads.txt must use /input/... paths: $container_path" >&2
                    exit 1
                    ;;
            esac
            host_path="$VS_INPUT_DIR/${container_path#/input/}"
            [[ -r "$host_path" ]] || {
                echo "ERROR: input is missing: $host_path" >&2
                exit 1
            }
        done < "$VS_READS_FILE"

        docker image inspect "$VS_IMAGE" >/dev/null 2>&1 || {
            echo "ERROR: image not found: $VS_IMAGE" >&2
            exit 1
        }

        docker rm -f "$VS_CONTAINER" >/dev/null 2>&1 || true
        rm -f "$VS_STATE_DIR/FAILED"

        docker run -d \
            --name "$VS_CONTAINER" \
            --platform linux/amd64 \
            --cpus="$VS_CPUS" \
            --memory="$VS_MEMORY" \
            --shm-size=4g \
            --user "$(id -u):$(id -g)" \
            --entrypoint /bin/bash \
            -e HOME=/tmp \
            -e VIRUSSEEKER_HOME=/opt/virusseeker \
            -e VS_DB_PROFILE="$VS_DB_PROFILE" \
            -e VS_CONFIG_FILE=/opt/virusseeker/scripts/VS.config \
            -e VS_CFG_FILE=/opt/virusseeker/scripts/VS.cfg \
            -e VS_LOCAL_ARRAY_JOBS=1 \
            -e VS_SBATCH_STATE_DIR=/state \
            -e FAMDB_DATA_DIR=/database/famdb \
            -e VS_REPEATMASKER_LIBRARY="${VS_REPEATMASKER_LIBRARY:-}" \
            -v "$VS_INPUT_DIR:/input:ro" \
            -v "$VS_OUTPUT_DIR:/opt/virusseeker/VS_output" \
            -v "$VS_DATABASE_DIR:/database:ro" \
            -v "$VS_SCRATCH_DIR:/scratch" \
            -v "$VS_STATE_DIR:/state" \
            --mount "type=bind,src=$VS_READS_FILE,dst=/work/reads.txt,readonly" \
            --mount "type=bind,src=$VS_CONFIG_FILE,dst=/opt/virusseeker/scripts/VS.config,readonly" \
            --mount "type=bind,src=$VS_CONFIG_FILE,dst=/opt/virusseeker/scripts/VS_Read_Counter/VS.config,readonly" \
            --mount "type=bind,src=$VS_CFG_FILE,dst=/opt/virusseeker/scripts/VS.cfg,readonly" \
            "$VS_IMAGE" \
            -c "
                set -Eeuo pipefail

                # Verify the exact flat VS.config seen by the current
                # parse_blast_Virus_noBP.py before starting queue.sh.
                /opt/conda/envs/vs/bin/python - <<'PY_CONTAINER_CFG'
from pathlib import Path

config_path = Path('/opt/virusseeker/scripts/VS.config')
flat = {}
with config_path.open('r', encoding='utf-8') as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line.startswith('#') or '=' not in line:
            continue
        key, value = line.split('=', 1)
        flat[key.strip()] = value.strip()

for key in ('vhunter', 'ncbi_taxdb', 'virus_family_genomes'):
    if not flat.get(key, ''):
        raise SystemExit(f'ERROR: mounted {config_path} missing flat key {key}')
    if not Path(flat[key]).exists():
        raise SystemExit(f'ERROR: configured path for {key} does not exist: {flat[key]}')

print('Mounted VS.config parser preflight OK:', config_path)
print('vhunter =', flat['vhunter'])
print('ncbi_taxdb =', flat['ncbi_taxdb'])
print('virus_family_genomes =', flat['virus_family_genomes'])
PY_CONTAINER_CFG

                grep -q 'bbsplit.sh path=.*bbmap_viral_ref_index' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing stage-24 BBSplitter index patch' >&2
                    exit 1
                }

                # Validate the positive stage-24 timestamp patches instead of
                # grepping globally for TIMEFILE. runVS.pl contains unrelated
                # TIMEFILE references elsewhere, so a global negative grep is
                # a false-positive for this runtime preflight.
                grep -q 'j24_map_to_viral_ref.time.txt' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing j24 timestamp patch' >&2
                    exit 1
                }
                grep -q 'j24b_map_to_viral_ref.time.txt' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing j24b timestamp patch' >&2
                    exit 1
                }

                grep -q 'bbmap.sh path=.*bbmap_host_index' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing host-removal BBMap index patch' >&2
                    exit 1
                }

                # Do not match literal shell-variable syntax for IN1/IN2 here.
                # This code is embedded in a double-quoted docker command, so
                # unescaped dollar-variable text would be expanded by the
                # outer shell before Docker starts.
                grep -E -q 'if \[ -s .*IN1.*\] && \[ -s .*IN2.*\]' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing metaSPAdes IN1/IN2 guard patch' >&2
                    exit 1
                }

                grep -q 'No paired reads remain after assembly; skipping PEAR cleanly.' /opt/virusseeker/scripts/runVS.pl || {
                    echo 'ERROR: mounted runVS.pl missing stage-6 empty-input PEAR patch' >&2
                    exit 1
                }

                echo 'Mounted runVS.pl stage-24 patch: OK'
                echo 'Mounted runVS.pl host-removal + assembly patches: OK'
                echo 'Mounted runVS.pl stage-6 empty-input PEAR patch: OK'

                test -s /opt/virusseeker/scripts/VS_Read_Counter/VS.config || {
                    echo 'ERROR: VS_Read_Counter/VS.config is not mounted' >&2
                    exit 1
                }
                cmp -s /opt/virusseeker/scripts/VS.config /opt/virusseeker/scripts/VS_Read_Counter/VS.config || {
                    echo 'ERROR: VS_Read_Counter/VS.config differs from scripts/VS.config' >&2
                    exit 1
                }
                echo 'Mounted VS_Read_Counter/VS.config: OK'

                cd /opt/virusseeker
                exec bash queue.sh \
                  -r /work/reads.txt \
                  -o \"$VS_RUN_NAME\" \
                  -d \"$VS_HOST_DB\" \
                  -a \"$VS_READ_TYPE\" \
                  -k \"$VS_ANALYSIS_MODE\" \
                  -m \"$VS_MEMORY\" \
                  -t \"$VS_CPUS\" \
                  ${VS_EXTRA_FLAGS:-}
            "

        echo "Started Docker container: $VS_CONTAINER"
        echo "Run directory: $RUN_DIR"
        echo "Monitor: docker logs -f $VS_CONTAINER"
        ;;

    slurm)
        command -v sbatch >/dev/null 2>&1 || {
            echo "ERROR: sbatch is not available." >&2
            exit 1
        }
        [[ -d "$VS_NATIVE_VS_HOME" ]] || {
            echo "ERROR: missing VS_NATIVE_VS_HOME: $VS_NATIVE_VS_HOME" >&2
            exit 1
        }

        # Native reads.txt contains real host paths.
        while IFS= read -r host_path; do
            host_path="${host_path%$'\r'}"
            [[ -z "$host_path" ]] && continue
            [[ -r "$host_path" ]] || {
                echo "ERROR: input is missing: $host_path" >&2
                exit 1
            }
        done < "$VS_READS_FILE"

        # The native scripts expect these fixed filenames in scripts/.
        for cfg in VS.config VS.cfg; do
            src="$VS_CONFIG_DIR/$cfg"
            dst="$VS_NATIVE_VS_HOME/scripts/$cfg"
            if [[ -e "$dst" && ! -L "$dst" && ! -e "$dst.repo-default" ]]; then
                cp -p "$dst" "$dst.repo-default"
            fi
            ln -sfn "$src" "$dst"
        done

        export VIRUSSEEKER_HOME="$VS_NATIVE_VS_HOME"
        export VS_CONFIG_FILE="$VS_NATIVE_VS_HOME/scripts/VS.config"
        export VS_CFG_FILE="$VS_NATIVE_VS_HOME/scripts/VS.cfg"
        export VS_DB_PROFILE="$VS_DB_PROFILE"

        mkdir -p "$VS_NATIVE_VS_HOME/scripts/VS_Read_Counter"
        ln -sfn "$VS_CONFIG_FILE" "$VS_NATIVE_VS_HOME/scripts/VS_Read_Counter/VS.config"
        export FAMDB_DATA_DIR="$VS_DATABASE_DIR/famdb"
        export VS_REPEATMASKER_LIBRARY="${VS_REPEATMASKER_LIBRARY:-}"

        # Validate the exact native VS.config used by the current Python parser.
        "$VS_NATIVE_CONDA_PREFIX/envs/vs/bin/python" - "$VS_CONFIG_FILE" <<'PY_NATIVE_CFG'
import sys

config_path = sys.argv[1]
flat = {}
with open(config_path, "r", encoding="utf-8") as fh:
    for raw in fh:
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        flat[key.strip()] = value.strip()

for key in ("vhunter", "ncbi_taxdb", "virus_family_genomes"):
    if not flat.get(key, ""):
        raise SystemExit(f"ERROR: native VS.config missing flat key {key}")

print("Native VS.config parser preflight OK:", config_path)
print("vhunter =", flat["vhunter"])
print("ncbi_taxdb =", flat["ncbi_taxdb"])
print("virus_family_genomes =", flat["virus_family_genomes"])
PY_NATIVE_CFG

        cd "$VS_NATIVE_VS_HOME"
        exec bash queue.sh \
          -r "$VS_READS_FILE" \
          -o "$VS_RUN_NAME" \
          -d "$VS_HOST_DB" \
          -a "$VS_READ_TYPE" \
          -k "$VS_ANALYSIS_MODE" \
          -m "$VS_MEMORY" \
          -t "$VS_CPUS" \
          ${VS_EXTRA_FLAGS:-}
        ;;

    *)
        echo "ERROR: unsupported VS_EXECUTION_MODE=$VS_EXECUTION_MODE" >&2
        exit 1
        ;;
esac
