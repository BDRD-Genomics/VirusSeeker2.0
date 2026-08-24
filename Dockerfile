# syntax=docker/dockerfile:1.7

ARG BASE_IMAGE=continuumio/miniconda3:24.11.1-0

FROM ${BASE_IMAGE} AS common

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
USER root

ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    CONDA_DIR=/opt/conda \
    CONDA_SOLVER=libmamba \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1
RUN --mount=type=cache,id=vs-apt-lists,target=/var/lib/apt/lists,sharing=locked \
    --mount=type=cache,id=vs-apt-cache,target=/var/cache/apt,sharing=locked \
    apt-get update \
    && apt-get install -y --no-install-recommends \
      bash \
      bzip2 \
      ca-certificates \
      coreutils \
      curl \
      findutils \
      gawk \
      git \
      grep \
      gzip \
      make \
      procps \
      sed \
      tar \
      unzip \
      wget \
      xz-utils
RUN cat > /opt/conda/.condarc <<'EOF_CONDARC'
channels:
  - conda-forge
  - bioconda
channel_priority: strict
solver: libmamba
auto_update_conda: false
show_channel_urls: true
EOF_CONDARC

RUN conda config --show solver \
    && conda config --show channel_priority \
    && conda --version

FROM common AS vs-env-builder

COPY docker/environment-vs.yml /tmp/environment-vs.yml

# perl-switch is no longer available in the selected channels. The environment
# includes perl-app-cpanminus; Switch 2.17 is installed from CPAN below.
RUN sed -i -E \
      '/^[[:space:]]*-[[:space:]]*perl-switch([[:space:]]*(#.*)?)?$/d' \
      /tmp/environment-vs.yml \
    && ! grep -nE \
      '^[[:space:]]*-[[:space:]]*perl-switch([[:space:]]|$)' \
      /tmp/environment-vs.yml

RUN --mount=type=cache,id=vs-conda-pkgs,target=/opt/conda/pkgs,sharing=locked \
    conda env create \
      --solver=libmamba \
      --file /tmp/environment-vs.yml \
    && test -x /opt/conda/envs/vs/bin/perl \
    && test -x /opt/conda/envs/vs/bin/cpanm

RUN --mount=type=cache,id=vs-cpan-cache,target=/root/.cpanm,sharing=locked \
    PATH="/opt/conda/envs/vs/bin:${PATH}" \
    PERL_MM_USE_DEFAULT=1 \
    /opt/conda/envs/vs/bin/cpanm --notest Switch@2.17 \
    && /opt/conda/envs/vs/bin/perl \
      -MSwitch \
      -e 'print "Switch.pm OK\n"'


FROM common AS dragonflye-env-builder

COPY docker/environment-dragonflye.yml /tmp/environment-dragonflye.yml

RUN --mount=type=cache,id=dragonflye-conda-pkgs,target=/opt/conda/pkgs,sharing=locked \
    conda env create \
      --solver=libmamba \
      --file /tmp/environment-dragonflye.yml \
    && test -x /opt/conda/envs/dragonflye/bin/dragonflye

FROM common AS repeatmasker-env-builder

COPY docker/environment-repeatmasker.yml /tmp/environment-repeatmasker.yml

RUN python - <<'PY_REPEATMASKER_YAML'
from pathlib import Path
import re

path = Path('/tmp/environment-repeatmasker.yml')
text = path.read_text()

if not re.search(r'^\s*-\s*h5py(?:\s|$)', text, flags=re.MULTILINE):
    lines = text.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if re.match(r'^\s*dependencies:\s*$', line):
            indent = re.match(r'^(\s*)', line).group(1) + '  '
            lines.insert(index + 1, f'{indent}- h5py\n')
            break
    else:
        raise SystemExit('dependencies: section not found in environment-repeatmasker.yml')
    path.write_text(''.join(lines))
PY_REPEATMASKER_YAML

RUN --mount=type=cache,id=repeatmasker-conda-pkgs,target=/opt/conda/pkgs,sharing=locked \
    conda env create \
      --solver=libmamba \
      --file /tmp/environment-repeatmasker.yml \
    && test -x /opt/conda/envs/repeatmasker/bin/RepeatMasker \
    && /opt/conda/envs/repeatmasker/bin/python -c \
      'import h5py; print("RepeatMasker h5py:", h5py.__version__)'

FROM common AS runtime

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
USER root

ARG VS_GIT_REF=unknown
ARG BUILD_DATE=unknown

ENV VIRUSSEEKER_HOME=/opt/virusseeker \
    VS_SBATCH_STATE_DIR=/tmp/vs-sbatch-state \
    VS_SBATCH_ARRAY_JOBS=1 \
    TMPDIR=/tmp

COPY --from=vs-env-builder /opt/conda/envs/vs /opt/conda/envs/vs
COPY --from=dragonflye-env-builder /opt/conda/envs/dragonflye /opt/conda/envs/dragonflye
COPY --from=repeatmasker-env-builder /opt/conda/envs/repeatmasker /opt/conda/envs/repeatmasker

RUN FAMDB_SCRIPT="$(find /opt/conda/envs/repeatmasker/share -type f -name famdb.py | head -n 1)" \
    && test -n "${FAMDB_SCRIPT}" \
    && FAMDB_INSTALL_DIR="$(dirname "${FAMDB_SCRIPT}")" \
    && sed -i '1c#!/opt/conda/envs/repeatmasker/bin/python' "${FAMDB_SCRIPT}" \
    && printf '%s\n' \
      '[famdb]' \
      'FAMDB_DATA_DIR = /database/famdb' \
      > "${FAMDB_INSTALL_DIR}/famdb.conf" \
    && /opt/conda/envs/repeatmasker/bin/python -c \
      'import h5py; print("RepeatMasker runtime h5py:", h5py.__version__)'
COPY . /opt/virusseeker

RUN chmod 0755 /opt/virusseeker/queue.sh \
    && bash -n /opt/virusseeker/queue.sh \
    && grep -q 'VIRUSSEEKER_HOME' /opt/virusseeker/queue.sh \
    && ! grep -qE \
         ';[[:space:]]*;;[[:space:]]*#.*removes[[:space:]]+human' \
         /opt/virusseeker/queue.sh \
    && grep -qE 'cp[[:space:]]+-fL[[:space:]].*R1' /opt/virusseeker/queue.sh \
    && grep -qE 'cp[[:space:]]+-fL[[:space:]].*R2' /opt/virusseeker/queue.sh \
    && grep -qE 'cp[[:space:]]+-fL[[:space:]].*\.fastq\.gz' /opt/virusseeker/queue.sh

RUN ! grep -q '\.unmapped_pe\.pear\.stitched\.fastq\.gz' \
      /opt/virusseeker/scripts/runVS.pl \
    && ! grep -q '\.unmapped\.pe\.pear\.stitched\.fastq\.gz' \
      /opt/virusseeker/scripts/runVS.pl \
    && /opt/conda/envs/vs/bin/perl \
      -c /opt/virusseeker/scripts/runVS.pl

COPY --chmod=0755 docker/sbatch-local.py /usr/local/bin/sbatch
COPY --chmod=0755 docker/test-sbatch-local.sh /usr/local/bin/test-sbatch-local
COPY docker/VS.config.container.example \
     /opt/virusseeker/scripts/VS.config.container.example


RUN mkdir -p \
      /path/to/VS \
      /data \
      /opt/virusseeker/VS_output \
      /opt/virusseeker/third_party/RepeatMasker \
      /opt/virusseeker/manifest \
      /database/famdb \
      /scratch/mmseqs \
      /scratch/diamond \
      /work \
      /tmp/vs-sbatch-state \
    && ln -sfn /opt/virusseeker /path/to/VS/dir \
    && ln -sfn /opt/virusseeker /data/VirusSeeker2.0 \
    && chmod 0777 \
      /opt/virusseeker/VS_output \
      /scratch/mmseqs \
      /scratch/diamond \
      /work \
      /tmp/vs-sbatch-state


RUN printf '%s\n' \
      '#!/usr/bin/env bash' \
      'set -euo pipefail' \
      'REPEATMASKER_ENV=/opt/conda/envs/repeatmasker' \
      'export PATH="${REPEATMASKER_ENV}/bin:${PATH}"' \
      'export PYTHONNOUSERSITE=1' \
      'export FAMDB_DATA_DIR="${FAMDB_DATA_DIR:-/database/famdb}"' \
      'if [[ -n "${VS_REPEATMASKER_LIBRARY:-}" ]]; then' \
      '  RM_LIB="${VS_REPEATMASKER_LIBRARY}"' \
      'elif [[ -s /database/test/repeat_library.fa ]]; then' \
      '  RM_LIB="/database/test/repeat_library.fa"' \
      'else' \
      '  RM_LIB=""' \
      'fi' \
      'if [[ -n "${RM_LIB}" ]]; then' \
      '  echo "RepeatMasker custom library: ${RM_LIB}" >&2' \
      '  export FAMDB_DIR=""' \
      '  unset FAMDB_DATA_DIR || true' \
      '  exec "${REPEATMASKER_ENV}/bin/RepeatMasker" -lib "${RM_LIB}" "$@"' \
      'else' \
      '  exec "${REPEATMASKER_ENV}/bin/RepeatMasker" "$@"' \
      'fi' \
      > /opt/virusseeker/third_party/RepeatMasker/RepeatMasker \
    && chmod 0755 \
      /opt/virusseeker/third_party/RepeatMasker/RepeatMasker \
    && find /opt/virusseeker/scripts /opt/virusseeker/bin \
      -type f \( -name '*.sh' -o -name '*.pl' -o -name '*.py' \) \
      -exec chmod a+rx {} +

ENV PATH="/usr/local/bin:/opt/conda/envs/vs/bin:/opt/conda/envs/dragonflye/bin:/opt/conda/envs/repeatmasker/bin:/opt/conda/bin:/opt/virusseeker/scripts:/opt/virusseeker/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin" \
    CONDA_DEFAULT_ENV=vs \
    CONDA_PREFIX=/opt/conda/envs/vs


RUN cat > /etc/profile.d/virusseeker-path.sh <<'EOF_VS_PATH'
export PATH="/usr/local/bin:/opt/conda/envs/vs/bin:/opt/conda/envs/dragonflye/bin:/opt/conda/envs/repeatmasker/bin:/opt/conda/bin:/opt/virusseeker/scripts:/opt/virusseeker/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin"
export CONDA_DEFAULT_ENV="vs"
export CONDA_PREFIX="/opt/conda/envs/vs"
EOF_VS_PATH

RUN chmod 0644 /etc/profile.d/virusseeker-path.sh

RUN cp \
      /opt/virusseeker/scripts/runVS.pl \
      /opt/virusseeker/scripts/runVS.pl.original \
    && cmp \
      /opt/virusseeker/scripts/runVS.pl \
      /opt/virusseeker/scripts/runVS.pl.original \
    && /opt/conda/envs/vs/bin/perl \
      -c /opt/virusseeker/scripts/runVS.pl \
    && /usr/local/bin/sbatch --version \
    && test "$(readlink -f /path/to/VS/dir)" = "/opt/virusseeker" \
    && test "$(readlink -f /data/VirusSeeker2.0)" = "/opt/virusseeker" \
    && test -r /opt/virusseeker/queue.sh \
    && bash -n /opt/virusseeker/queue.sh

RUN conda list --prefix /opt/conda/envs/vs --explicit \
      > /opt/virusseeker/manifest/vs-explicit.txt \
    && conda list --prefix /opt/conda/envs/dragonflye --explicit \
      > /opt/virusseeker/manifest/dragonflye-explicit.txt \
    && conda list --prefix /opt/conda/envs/repeatmasker --explicit \
      > /opt/virusseeker/manifest/repeatmasker-explicit.txt \
    && printf '%s\n' \
      'Switch 2.17 installed from CPAN with cpanm' \
      > /opt/virusseeker/manifest/cpan-modules.txt

RUN cat > /usr/local/bin/virusseeker-entrypoint <<'EOF_ENTRYPOINT'
#!/usr/bin/env bash
set -euo pipefail

case "${1:-bash}" in
  runVS)
    shift
    cd /opt/virusseeker/scripts
    exec /opt/conda/envs/vs/bin/perl runVS.pl "$@"
    ;;
  queue)
    shift
    cd /opt/virusseeker
    exec bash queue.sh "$@"
    ;;
  test-sbatch)
    shift
    exec /usr/local/bin/test-sbatch-local "$@"
    ;;
  bash|sh)
    shell="$1"
    shift
    exec "/bin/${shell}" "$@"
    ;;
  *)
    exec "$@"
    ;;
esac
EOF_ENTRYPOINT

RUN chmod 0755 /usr/local/bin/virusseeker-entrypoint

RUN cat > /usr/local/bin/virusseeker-healthcheck <<'EOF_HEALTHCHECK'
#!/usr/bin/env bash
set -euo pipefail

command -v sbatch >/dev/null
test -x /opt/conda/envs/vs/bin/perl
command -v python >/dev/null
command -v mmseqs >/dev/null
command -v diamond >/dev/null
command -v bbduk.sh >/dev/null

test -x /opt/conda/envs/dragonflye/bin/dragonflye
test -x /opt/conda/envs/repeatmasker/bin/RepeatMasker
test -x /opt/virusseeker/third_party/RepeatMasker/RepeatMasker
/opt/conda/envs/repeatmasker/bin/python -c 'import h5py'
test -f "$(find /opt/conda/envs/repeatmasker/share -type f -name famdb.conf | head -n 1)"

cmp \
  /opt/virusseeker/scripts/runVS.pl \
  /opt/virusseeker/scripts/runVS.pl.original

/opt/conda/envs/vs/bin/perl -MConfig::Simple -e '1'
/opt/conda/envs/vs/bin/perl -MXML::Simple -e '1'
/opt/conda/envs/vs/bin/perl -MSwitch -e '1'
/opt/conda/envs/vs/bin/perl -c /opt/virusseeker/scripts/runVS.pl >/dev/null
sbatch --version >/dev/null
EOF_HEALTHCHECK

RUN chmod 0755 /usr/local/bin/virusseeker-healthcheck


RUN find /opt/virusseeker \
        -type d \
        -exec chmod 0755 {} + \
    && find /opt/virusseeker/scripts \
        -type f \
        -name '*.pl' \
        -exec chmod 0755 {} + \
    && find /opt/virusseeker \
        -maxdepth 1 \
        -type f \
        -name '*.sh' \
        -exec chmod 0755 {} +


LABEL org.opencontainers.image.title="VirusSeeker 2.0 local sbatch runtime" \
      org.opencontainers.image.description="VirusSeeker 2.0 Docker-compatible pipeline with a single-node sbatch compatibility executor" \
      org.opencontainers.image.source="https://github.com/BDRD-Genomics/VirusSeeker2.0" \
      org.opencontainers.image.revision="${VS_GIT_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

RUN printf '%s\n' \
      "VirusSeeker source revision: ${VS_GIT_REF}" \
      "Build date: ${BUILD_DATE}" \
      > /opt/virusseeker/manifest/build.txt

VOLUME ["/database", "/scratch", "/opt/virusseeker/VS_output"]
WORKDIR /work

HEALTHCHECK --interval=30s --timeout=20s --start-period=10s --retries=3 \
  CMD ["/usr/local/bin/virusseeker-healthcheck"]

ENTRYPOINT ["/usr/local/bin/virusseeker-entrypoint"]
CMD ["bash"]