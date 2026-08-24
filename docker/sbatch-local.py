#!/opt/conda/envs/vs/bin/python
"""Synchronous sbatch compatibility shim for VirusSeeker 2.0.

"""

from __future__ import annotations

import fcntl
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

VERSION = "slurm 23.11.0 (VirusSeeker local sbatch compatibility)"
DIRECTIVE_RE = re.compile(r"^\s*#SBATCH\s+--([A-Za-z0-9_-]+)(?:[=\s]+(.*?))?\s*$")


MEMORY_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([KMGTkmgt]?)\s*(?:[Bb])?\s*$")


def memory_to_mb(value: Optional[str]) -> str:
    """Convert a Slurm-style memory value to integer megabytes.
    """
    if value is None:
        return "0"

    raw = str(value).strip()
    if not raw:
        return "0"

    match = MEMORY_RE.fullmatch(raw)
    if not match:
        die(f"unsupported Slurm memory value: {value}")

    number = float(match.group(1))
    unit = match.group(2).upper()

    factors = {
        "": 1,
        "M": 1,
        "K": 1 / 1024,
        "G": 1024,
        "T": 1024 * 1024,
    }

    megabytes = int(-(-number * factors[unit] // 1))
    return str(max(1, megabytes))


def die(message: str, code: int = 1) -> "NoReturn":
    print(f"sbatch-local: {message}", file=sys.stderr)
    raise SystemExit(code)


def parse_cli(argv: List[str]) -> Tuple[Path, List[str]]:
    if not argv:
        die("missing batch script")
    if argv[0] in {"--version", "-V"}:
        print(VERSION)
        raise SystemExit(0)

    # VirusSeeker calls: sbatch /path/to/script.sh
    # Accept a few harmless common flags so the shim remains easy to inspect/use.
    script: Optional[Path] = None
    script_args: List[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if script is not None:
            script_args.append(arg)
        elif arg in {"--wait", "--parsable"}:
            pass
        elif arg.startswith("-"):
            # Current VirusSeeker does not pass sbatch CLI options. Fail loudly
            # rather than silently misinterpreting an unsupported invocation.
            die(f"unsupported command-line option: {arg}")
        else:
            script = Path(arg)
        i += 1

    if script is None:
        die("missing batch script")
    if not script.is_file():
        die(f"batch script does not exist: {script}")
    return script.resolve(), script_args


def parse_directives(script: Path) -> Dict[str, str]:
    directives: Dict[str, str] = {}
    with script.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = DIRECTIVE_RE.match(line)
            if match:
                key = match.group(1)
                value = (match.group(2) or "").strip().strip('"').strip("'")
                directives[key] = value
    return directives


def expand_array(spec: Optional[str]) -> List[Optional[int]]:
    if not spec:
        return [None]

    # Ignore SLURM's optional %concurrency suffix; local execution is controlled
    # by VS_SBATCH_ARRAY_JOBS and is sequential by default.
    raw = spec.split("%", 1)[0]
    tasks: List[int] = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        match = re.fullmatch(r"(-?\d+)-(-?\d+)(?::(\d+))?", part)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            step = int(match.group(3) or "1")
            if step <= 0:
                die(f"invalid array step in: {spec}")
            if end < start:
                die(f"descending arrays are unsupported: {spec}")
            tasks.extend(range(start, end + 1, step))
        elif re.fullmatch(r"-?\d+", part):
            tasks.append(int(part))
        else:
            die(f"unsupported array specification: {spec}")

    if not tasks:
        die(f"empty array specification: {spec}")

    limit_raw = os.environ.get("VS_SBATCH_ARRAY_LIMIT", "0")
    try:
        limit = int(limit_raw)
    except ValueError:
        die(f"VS_SBATCH_ARRAY_LIMIT must be an integer, got: {limit_raw}")
    if limit > 0 and len(tasks) > limit:
        print(
            f"sbatch-local: smoke-test cap enabled; executing first {limit} "
            f"of {len(tasks)} array tasks from {spec}",
            file=sys.stderr,
        )
        tasks = tasks[:limit]
    return [int(task) for task in tasks]


def allocate_job_id(state_dir: Path) -> int:
    state_dir.mkdir(parents=True, exist_ok=True)
    counter = state_dir / "next_job_id"
    with counter.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        raw = handle.read().strip()
        current = int(raw) if raw else 1000
        job_id = current + 1
        handle.seek(0)
        handle.truncate()
        handle.write(str(job_id))
        handle.flush()
        os.fsync(handle.fileno())
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return job_id


def render_path(template: str, job_id: int, task_id: Optional[int], job_name: str) -> Path:
    task_token = str(task_id if task_id is not None else 0)
    rendered = template
    rendered = rendered.replace("%%", "\0")
    rendered = rendered.replace("%A", str(job_id))
    rendered = rendered.replace("%a", task_token)
    rendered = rendered.replace("%j", str(job_id))
    rendered = rendered.replace("%x", job_name)
    rendered = rendered.replace("\0", "%")
    return Path(rendered)


def task_environment(
    directives: Dict[str, str], job_id: int, task_id: Optional[int], submit_dir: Path
) -> Dict[str, str]:
    env = os.environ.copy()
    job_name = directives.get("job-name") or f"job{job_id}"
    env.update(
        {
            "SLURM_JOB_ID": str(job_id),
            "SLURM_JOB_NAME": job_name,
            "SLURM_SUBMIT_DIR": str(submit_dir),
            "SLURM_CPUS_PER_TASK": directives.get("cpus-per-task", "1"),
            "SLURM_MEM_PER_NODE": memory_to_mb(directives.get("mem", "0")),
            "SLURM_JOB_PARTITION": directives.get("partition", "local"),
        }
    )
    if task_id is not None:
        env["SLURM_ARRAY_JOB_ID"] = str(job_id)
        env["SLURM_ARRAY_TASK_ID"] = str(task_id)
    else:
        env.pop("SLURM_ARRAY_JOB_ID", None)
        env.pop("SLURM_ARRAY_TASK_ID", None)
    return env


def open_log(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="utf-8")


def run_task(
    script: Path,
    script_args: List[str],
    directives: Dict[str, str],
    job_id: int,
    task_id: Optional[int],
    submit_dir: Path,
) -> int:
    job_name = directives.get("job-name") or f"job{job_id}"
    default_out = "slurm-%A_%a.out" if task_id is not None else "slurm-%j.out"
    out_template = directives.get("output") or default_out
    err_template = directives.get("error") or out_template
    out_path = render_path(out_template, job_id, task_id, job_name)
    err_path = render_path(err_template, job_id, task_id, job_name)

    env = task_environment(directives, job_id, task_id, submit_dir)
    cmd = ["/bin/bash", str(script), *script_args]

    # If stdout and stderr target the same path, use one file descriptor.
    if out_path == err_path:
        with open_log(out_path) as combined:
            completed = subprocess.run(
                cmd,
                cwd=submit_dir,
                env=env,
                stdout=combined,
                stderr=subprocess.STDOUT,
                check=False,
            )
    else:
        with open_log(out_path) as stdout_handle, open_log(err_path) as stderr_handle:
            completed = subprocess.run(
                cmd,
                cwd=submit_dir,
                env=env,
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
            )
    return completed.returncode


def main(argv: List[str]) -> int:
    script, script_args = parse_cli(argv)
    directives = parse_directives(script)
    state_dir = Path(os.environ.get("VS_SBATCH_STATE_DIR", "/tmp/vs-sbatch-state"))
    failed_sentinel = state_dir / "FAILED"

    if failed_sentinel.exists():
        die(
            f"a previous local job failed; remove {failed_sentinel} after reviewing logs",
            code=1,
        )

    job_id = allocate_job_id(state_dir)
    submit_dir = Path.cwd().resolve()
    tasks = expand_array(directives.get("array"))

    for task_id in tasks:
        rc = run_task(
            script=script,
            script_args=script_args,
            directives=directives,
            job_id=job_id,
            task_id=task_id,
            submit_dir=submit_dir,
        )
        if rc != 0:
            failed_sentinel.parent.mkdir(parents=True, exist_ok=True)
            failed_sentinel.write_text(
                f"job_id={job_id}\nscript={script}\ntask_id={task_id}\nexit_code={rc}\n",
                encoding="utf-8",
            )
            die(
                f"job {job_id} task {task_id if task_id is not None else 'N/A'} "
                f"failed with exit code {rc}; see SBATCH output/error logs",
                code=rc,
            )

    print(f"Submitted batch job {job_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
