# noxfile.py
"""
Run the model's bash entrypoint in a controlled session (no interactive activation).
Usage examples:

  # Using absolute env prefix (recommended)
  nox -s run_model -- \
    --env-prefix /usr/share/miniconda/envs/eos9o72 \
    --run-dir /home/runner/eos/temp/eos9o72/model/framework \
    --input /home/runner/eos/temp/eos9o72/model/framework/examples/run_input.csv \
    --out /tmp/bash_output.csv

  # Or using env name (falls back to 'conda info --base' to resolve prefix)
  nox -s run_model -- \
    --env-name eos9o72 \
    --run-dir /home/runner/eos/temp/eos9o72/model/framework \
    --input /home/runner/eos/temp/eos9o72/model/framework/examples/run_input.csv \
    --out /tmp/bash_output.csv
"""

import os
import shlex
import subprocess
from pathlib import Path

import nox


def _get_arg(args, key, default=None):
    if key in args:
        i = args.index(key)
        return args[i + 1] if i + 1 < len(args) else default
    return default


def _conda_base():
    # Try CONDA_PREFIX or CONDA_EXE, else conda info --base
    conda_exe = os.environ.get("CONDA_EXE") or "conda"
    # If CONDA_EXE points to .../bin/conda, derive base
    if os.path.basename(conda_exe) == "conda":
        base_guess = Path(conda_exe).resolve().parents[1] if "miniconda" in conda_exe or "anaconda" in conda_exe else None
        if base_guess and base_guess.exists():
            return str(base_guess)
    try:
        out = subprocess.check_output([conda_exe, "info", "--base"], text=True).strip()
        return out
    except Exception:
        return None


def _resolve_env_prefix(env_prefix, env_name):
    if env_prefix:
        return env_prefix
    if env_name:
        base = _conda_base()
        if not base:
            raise RuntimeError("Could not resolve conda base to compute env prefix from env name.")
        return str(Path(base) / "envs" / env_name)
    raise RuntimeError("You must provide --env-prefix or --env-name.")


@nox.session
def run_model(session: nox.Session):
    """
    Executes model/framework/run.sh with (.) INPUT OUTPUT inside a hermetic env.

    Priority:
      1) micromamba run -p <prefix>
      2) conda run      -p <prefix>
    """
    args = list(session.posargs)

    env_prefix = _get_arg(args, "--env-prefix")
    env_name = _get_arg(args, "--env-name")
    run_dir = _get_arg(args, "--run-dir")
    inp = _get_arg(args, "--input")
    out = _get_arg(args, "--out")

    if not run_dir or not inp or not out:
        raise SystemExit("Missing required args. Provide --run-dir, --input, --out and either --env-prefix or --env-name.")

    env_prefix = _resolve_env_prefix(env_prefix, env_name)

    run_dir = str(Path(run_dir).resolve())
    inp = str(Path(inp).resolve())
    out = str(Path(out).resolve())

    # Ensure run.sh is executable
    session.run("bash", "-c", f"cd {shlex.quote(run_dir)} && chmod +x ./run.sh", external=True)

    # Prefer micromamba if present; otherwise conda.
    if shutil.which("micromamba"):
        cmd = (
            f'micromamba run -p {shlex.quote(env_prefix)} '
            f'bash -c {shlex.quote(f"cd {run_dir} && ./run.sh . {inp} {out}")}'
        )
    else:
        conda_exe = os.environ.get("CONDA_EXE", "conda")
        cmd = (
            f'{shlex.quote(conda_exe)} run -p {shlex.quote(env_prefix)} '
            f'bash -c {shlex.quote(f"cd {run_dir} && ./run.sh . {inp} {out}")}'
        )

    # Show which interpreter will run
    show_py = (
        f'python -V && which python && echo "CONDA_PREFIX=${{CONDA_PREFIX:-}}"; '
        f'cd {shlex.quote(run_dir)} && ls -lah'
    )
    if "micromamba run" in cmd:
        session.run("bash", "-lc", f'micromamba run -p {shlex.quote(env_prefix)} bash -c {shlex.quote(show_py)}', external=True)
    else:
        conda_exe = os.environ.get("CONDA_EXE", "conda")
        session.run("bash", "-lc", f'{shlex.quote(conda_exe)} run -p {shlex.quote(env_prefix)} bash -c {shlex.quote(show_py)}', external=True)

    # Run the model
    session.log(f"Executing: {cmd}")
    session.run("bash", "-lc", cmd, external=True)

    # Verify output exists
    session.run("bash", "-lc", f'test -f {shlex.quote(out)} || (echo "ERROR: output not found: {out}" >&2; exit 1)', external=True)


# Utilities needed at top level
import shutil  # after function defs to keep nox import first
