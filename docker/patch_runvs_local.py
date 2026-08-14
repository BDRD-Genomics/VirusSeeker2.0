#!/usr/bin/env python3
"""Add a checked local/slurm submission helper to VirusSeeker runVS.pl."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

MARKER = "VIRUSSEEKER_LOCAL_EXECUTOR_PATCH_V1"

HELPER = r'''
#####################################################################################
# VIRUSSEEKER_LOCAL_EXECUTOR_PATCH_V1
# Container execution defaults to a blocking local submitter. Because each submitted
# stage completes before submit_job() returns, the existing stage order supplies the
# same afterok behavior without a scheduler. Set VS_EXECUTOR=slurm to retain sbatch.
if (($ENV{'VS_EXECUTOR'} // 'slurm') eq 'local' && $use_checkpoint) {
    print STDERR "[VirusSeeker] checkpoint mode is unavailable locally; continuing with checkpointing disabled.\n";
    $use_checkpoint = 0;
}

sub submit_job {
    my ($script_path) = @_;
    my $executor = $ENV{'VS_EXECUTOR'} // 'slurm';
    my $local_submit = $ENV{'VS_LOCAL_SUBMIT'} // '/usr/local/bin/vs-local-submit';
    die "Missing job script path\n" unless defined $script_path && length $script_path;

    if ($executor eq 'local') {
        print STDERR "[VirusSeeker] local execution: $script_path\n";
        system($local_submit, $script_path);
        if ($? == -1) {
            die "Unable to launch local submitter '$local_submit': $!\n";
        }
        if ($? & 127) {
            die sprintf("Local job died from signal %d: %s\n", ($? & 127), $script_path);
        }
        my $exit_code = $? >> 8;
        die "Local job failed with exit code $exit_code: $script_path\n" if $exit_code != 0;
        return "local";
    }

    if ($executor eq 'slurm') {
        my $output = qx{sbatch "$script_path"};
        my $exit_code = $? >> 8;
        die "sbatch failed with exit code $exit_code: $script_path\n" if $exit_code != 0;
        if ($output =~ /(\d+)\s*$/) {
            return $1;
        }
        die "Could not parse SLURM job id from: $output\n";
    }

    die "Unsupported VS_EXECUTOR '$executor'; expected local or slurm\n";
}
#####################################################################################
'''


def patch(path: Path) -> int:
    text = path.read_text(encoding="utf-8")
    if MARKER in text:
        print(f"Already patched: {path}")
        return 0

    # Replace the original `sbatch ... | awk` assignments while preserving which
    # last_jobid/current_job_file variable each stage uses.
    pattern = re.compile(
        r"\$(last_jobid[12])\s*=\s*`sbatch\s+\$job_files_dir/\$(current_job_file[12])[^`]*`;"
    )

    def replacement(match: re.Match[str]) -> str:
        return (
            f'${match.group(1)} = '
            f'submit_job("$job_files_dir/${match.group(2)}");'
        )

    patched, count = pattern.subn(replacement, text)
    if count < 10:
        raise RuntimeError(
            f"Expected to replace at least 10 sbatch submissions, replaced {count}. "
            "The upstream runVS.pl layout may have changed."
        )

    anchor = "####################################################################################\n# run the whole pipeline"
    if anchor not in patched:
        raise RuntimeError("Could not find run-the-whole-pipeline insertion anchor")
    patched = patched.replace(anchor, HELPER + "\n" + anchor, 1)

    remaining = [
        line for line in patched.splitlines()
        if "`sbatch " in line and not line.lstrip().startswith("#")
    ]
    if remaining:
        raise RuntimeError(f"Unpatched sbatch calls remain: {remaining[:3]}")

    path.write_text(patched, encoding="utf-8")
    print(f"Patched {path}: replaced {count} sbatch submissions")
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    patch(args.path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
