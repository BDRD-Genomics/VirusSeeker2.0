#!/usr/bin/env bash
set -euo pipefail

run_vs() {
  local executor="$1"
  shift
  export VS_EXECUTOR="$executor"
  exec perl /opt/virusseeker/scripts/runVS.pl "$@"
}

case "${1:-bash}" in
  healthcheck)
    shift
    exec /usr/local/bin/virusseeker-healthcheck "$@"
    ;;
  vs)
    shift
    exec conda run --no-capture-output --name vs "$@"
    ;;
  dragonflye-env)
    shift
    exec conda run --no-capture-output --name dragonflye "$@"
    ;;
  repeatmasker-env)
    shift
    exec conda run --no-capture-output --name repeatmasker "$@"
    ;;
  runVS|runVS-local)
    shift
    run_vs local "$@"
    ;;
  runVS-slurm)
    shift
    run_vs slurm "$@"
    ;;
  local-submit)
    shift
    exec /usr/local/bin/vs-local-submit "$@"
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
