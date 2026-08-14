#!/usr/bin/env bash
set -euo pipefail

quick=false
if [[ "${1:-}" == "--quick" ]]; then
  quick=true
fi

required=(
  perl python fastp bbduk.sh reformat.sh blastn bwa samtools bamToFastq
  spades.py unicycler diamond mmseqs minimap2 seqkit seqtk pear pigz
  cutadapt taxonkit fastqc vs-local-submit
)

failed=0
for exe in "${required[@]}"; do
  if ! command -v "$exe" >/dev/null 2>&1; then
    printf 'MISSING\t%s\n' "$exe" >&2
    failed=1
  elif ! $quick; then
    printf 'FOUND\t%s\t%s\n' "$exe" "$(command -v "$exe")"
  fi
done

for path in \
  /opt/conda/envs/dragonflye/bin/dragonflye \
  /opt/conda/envs/repeatmasker/bin/RepeatMasker \
  /opt/virusseeker/scripts/runVS.pl \
  /opt/virusseeker/scripts/runVS.pl.slurm-original \
  /opt/virusseeker/scripts/VS.config.container.example; do
  if [[ ! -e "$path" ]]; then
    printf 'MISSING\t%s\n' "$path" >&2
    failed=1
  elif ! $quick; then
    printf 'FOUND\t%s\n' "$path"
  fi
done

perl -MConfig::Simple -e 'print "Config::Simple OK\n"' >/dev/null
perl -MXML::Simple -e 'print "XML::Simple OK\n"' >/dev/null
grep -q 'VIRUSSEEKER_LOCAL_EXECUTOR_PATCH_V1' /opt/virusseeker/scripts/runVS.pl

if ! $quick; then
  vs-local-submit --self-test >/dev/null
fi

if (( failed != 0 )); then
  exit 1
fi

if ! $quick; then
  printf '\nVirusSeeker runtime healthcheck passed.\n'
  printf 'Default executor:        local (no SLURM required)\n'
  printf 'Main environment:        /opt/conda/envs/vs\n'
  printf 'Dragonflye environment:  /opt/conda/envs/dragonflye\n'
  printf 'RepeatMasker environment:/opt/conda/envs/repeatmasker\n'
fi
