#!/usr/bin/env python3
"""Run a generated VirusSeeker SLURM job script locally.

This intentionally implements only the SBATCH directives used by VirusSeeker:
--array, --cpus-per-task, --mem, --output, --error, and --job-name.
Dependencies are unnecessary because this submitter is blocking: a stage does
not return until all of its local array tasks finish successfully.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

VERSION = "1.0.0"


@dataclass(frozen=True)
class Directives:
    array: str | None = None
    cpus: str = "1"
    memory: str = "1g"
    output: str | None = None
    error: str | None = None
    job_name: str = "virusseeker"


def parse_directives(script: Path) -> Directives:
    values: dict[str, str] = {}
    pattern = re.compile(r"^\s*#SBATCH\s+--([A-Za-z0-9_-]+)(?:=|\s+)(.*?)\s*$")
    with script.open("r", encoding="utf-8") as handle:
        for line in handle:
            match = pattern.match(line)
            if not match:
                continue
            values[match.group(1)] = match.group(2).strip().strip('"').strip("'")
    return Directives(
        array=values.get("array"),
        cpus=values.get("cpus-per-task", os.environ.get("VS_LOCAL_CPUS", "1")),
        memory=values.get("mem", os.environ.get("VS_LOCAL_MEMORY", "1g")),
        output=values.get("output"),
        error=values.get("error"),
        job_name=values.get("job-name", script.stem),
    )


def parse_array(spec: str | None) -> tuple[list[int | None], int | None]:
    if not spec:
        return [None], None
    raw, percent, limit_text = spec.partition("%")
    scheduler_limit = int(limit_text) if percent and limit_text.isdigit() else None
    tasks: list[int] = []
    for chunk in raw.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        match = re.fullmatch(r"(-?\d+)(?:-(-?\d+)(?::(\d+))?)?", chunk)
        if not match:
            raise ValueError(f"Unsupported array specification: {spec!r}")
        start = int(match.group(1))
        if match.group(2) is None:
            tasks.append(start)
            continue
        end = int(match.group(2))
        step = int(match.group(3) or "1")
        if step <= 0:
            raise ValueError(f"Array step must be positive: {spec!r}")
        stop = end + (1 if end >= start else -1)
        signed_step = step if end >= start else -step
        tasks.extend(range(start, stop, signed_step))
    if not tasks:
        raise ValueError(f"Array specification produced no tasks: {spec!r}")
    return tasks, scheduler_limit


def expand_log_path(template: str | None, job_id: str, task_id: int | None, job_name: str) -> Path | None:
    if not template:
        return None
    task_token = str(task_id) if task_id is not None else "0"
    expanded = (
        template.replace("%A", job_id)
        .replace("%a", task_token)
        .replace("%j", job_id)
        .replace("%x", job_name.replace("/", "_"))
    )
    return Path(expanded)


def run_task(
    script: Path,
    directives: Directives,
    job_id: str,
    task_id: int | None,
) -> tuple[int | None, int]:
    env = os.environ.copy()
    env.update(
        {
            "SLURM_JOB_ID": job_id,
            "SLURM_ARRAY_JOB_ID": job_id,
            "SLURM_ARRAY_TASK_ID": str(task_id if task_id is not None else 0),
            "SLURM_CPUS_PER_TASK": directives.cpus,
            "SLURM_MEM_PER_NODE": directives.memory,
            "SLURM_SUBMIT_DIR": os.getcwd(),
            "VS_EXECUTOR": "local",
        }
    )

    stdout_path = expand_log_path(directives.output, job_id, task_id, directives.job_name)
    stderr_path = expand_log_path(directives.error, job_id, task_id, directives.job_name)
    if stdout_path:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    if stderr_path:
        stderr_path.parent.mkdir(parents=True, exist_ok=True)

    stdout_handle = stdout_path.open("wb") if stdout_path else None
    if stderr_path:
        stderr_handle = stderr_path.open("wb")
    elif stdout_handle:
        stderr_handle = stdout_handle
    else:
        stderr_handle = None

    label = f"task {task_id}" if task_id is not None else "job"
    print(
        f"[vs-local] starting {label}: {script} (cpus={directives.cpus}, mem={directives.memory})",
        file=sys.stderr,
        flush=True,
    )
    try:
        completed = subprocess.run(
            ["/bin/bash", str(script)],
            env=env,
            cwd=os.getcwd(),
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )
        return task_id, completed.returncode
    finally:
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle and stderr_handle is not stdout_handle:
            stderr_handle.close()


def choose_workers(task_count: int, scheduler_limit: int | None) -> int:
    requested = int(os.environ.get("VS_LOCAL_ARRAY_JOBS", "1"))
    if requested < 1:
        raise ValueError("VS_LOCAL_ARRAY_JOBS must be at least 1")
    if scheduler_limit is not None:
        requested = min(requested, scheduler_limit)
    return max(1, min(requested, task_count))


def execute(script: Path) -> str:
    if not script.is_file():
        raise FileNotFoundError(f"Job script does not exist: {script}")
    directives = parse_directives(script)
    tasks, scheduler_limit = parse_array(directives.array)
    workers = choose_workers(len(tasks), scheduler_limit)
    job_id = str(time.time_ns())

    failures: list[tuple[int | None, int]] = []
    if workers == 1:
        for task_id in tasks:
            item = run_task(script, directives, job_id, task_id)
            if item[1] != 0:
                failures.append(item)
                break
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(run_task, script, directives, job_id, task_id) for task_id in tasks]
            for future in concurrent.futures.as_completed(futures):
                item = future.result()
                if item[1] != 0:
                    failures.append(item)

    if failures:
        details = ", ".join(
            f"task {task if task is not None else 'job'} exit={code}" for task, code in failures
        )
        raise RuntimeError(f"Local VirusSeeker stage failed: {details}; script={script}")

    # Keep the same stdout shape expected by the original awk parser.
    print(f"Submitted batch job {job_id}")
    return job_id


def self_test() -> None:
    import tempfile

    with tempfile.TemporaryDirectory(prefix="vs-local-submit-") as tmp:
        root = Path(tmp)
        script = root / "array.sh"
        script.write_text(
            "#!/bin/bash\n"
            "#SBATCH --array=1-3\n"
            "#SBATCH --cpus-per-task=2\n"
            "#SBATCH --mem=4g\n"
            f"#SBATCH --output={root}/task.%A_%a.out\n"
            "echo ${SLURM_ARRAY_TASK_ID}:${SLURM_CPUS_PER_TASK}:${SLURM_MEM_PER_NODE}\n",
            encoding="utf-8",
        )
        old_jobs = os.environ.get("VS_LOCAL_ARRAY_JOBS")
        os.environ["VS_LOCAL_ARRAY_JOBS"] = "2"
        try:
            execute(script)
        finally:
            if old_jobs is None:
                os.environ.pop("VS_LOCAL_ARRAY_JOBS", None)
            else:
                os.environ["VS_LOCAL_ARRAY_JOBS"] = old_jobs
        outputs = sorted(root.glob("task.*_*.out"))
        if len(outputs) != 3:
            raise RuntimeError(f"Self-test expected 3 logs, found {len(outputs)}")
        observed = sorted(path.read_text(encoding="utf-8").strip() for path in outputs)
        expected = ["1:2:4g", "2:2:4g", "3:2:4g"]
        if observed != expected:
            raise RuntimeError(f"Self-test mismatch: {observed!r} != {expected!r}")
    print("vs-local-submit self-test passed", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("script", nargs="?", type=Path)
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        if args.script is None:
            parser.error("a job script is required unless --self-test is used")
        execute(args.script.resolve())
        return 0
    except Exception as exc:  # concise fatal error for Perl caller and users
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
