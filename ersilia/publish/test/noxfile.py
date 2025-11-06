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
import os, shlex, shutil
from pathlib import Path
import nox


@nox.session
def run_model(session: nox.Session):
    """Run the model's run.sh inside the proper conda/mamba env and print full logs."""
    args = list(session.posargs)
    env_prefix = _get_arg(args, "--env-prefix")
    env_name   = _get_arg(args, "--env-name")
    run_dir    = Path(_get_arg(args, "--run-dir")).resolve()
    inp        = Path(_get_arg(args, "--input")).resolve()
    out_final  = Path(_get_arg(args, "--out")).resolve()

    if not run_dir or not inp or not out_final:
        raise SystemExit("Missing required args: --run-dir, --input, --out plus --env-prefix/--env-name.")

    env_prefix = _resolve_env_prefix(env_prefix, env_name)
    out_in_run = run_dir / "ersilia_bash_output.csv"

    session.log(f"=== MODEL RUN DEBUG INFO ===")
    session.log(f"run_dir: {run_dir}")
    session.log(f"env_prefix: {env_prefix}")
    session.log(f"input file: {inp}")
    session.log(f"output (expected inside run_dir): {out_in_run}")
    session.log(f"final output copy target: {out_final}")

    # Ensure script is executable
    session.run("bash", "-c", f"cd {shlex.quote(str(run_dir))} && chmod +x ./run.sh", external=True)

    # Print input file content for verification
    session.log("\n=== INPUT FILE CONTENT ===")
    session.run("bash", "-c", f"echo '--- {inp} ---'; cat {shlex.quote(str(inp))} || echo '(no input found)'", external=True)

    # List run_dir contents before running
    session.log("\n=== CONTENTS OF run_dir BEFORE EXECUTION ===")
    session.run("bash", "-c", f"ls -lah {shlex.quote(str(run_dir))}", external=True)

    # Construct command to run inside env
    run_cmd = f'cd {run_dir} && ./run.sh . {inp} {out_in_run}'
    show_env = (
        f'python -V && which python && '
        f'echo "CONDA_PREFIX=${{CONDA_PREFIX:-}}" && '
        f'echo "--- run_dir: {run_dir} ---" && ls -lah {run_dir}'
    )

    def run_in_env(command: str):
        if shutil.which("micromamba"):
            session.run(
                "bash",
                "-lc",
                f'micromamba run -p {shlex.quote(env_prefix)} bash -c {shlex.quote(command)}',
                external=True,
            )
        else:
            conda_exe = os.environ.get("CONDA_EXE", "conda")
            session.run(
                "bash",
                "-lc",
                f'{shlex.quote(conda_exe)} run -p {shlex.quote(env_prefix)} bash -c {shlex.quote(command)}',
                external=True,
            )

    session.log("\n=== PYTHON ENVIRONMENT INSIDE CONDA ===")
    run_in_env(show_env)

    session.log("\n=== RUN.SH EXECUTION ===")
    run_in_env(run_cmd)

    # Show directory contents after run
    session.log("\n=== CONTENTS OF run_dir AFTER EXECUTION ===")
    session.run("bash", "-c", f"ls -lah {shlex.quote(str(run_dir))}", external=True)

    # If expected output missing, try to locate another CSV
    if not out_in_run.exists():
        session.log("\n⚠️  Expected output not found, searching for alternative CSVs...")
        candidates = list(sorted(run_dir.glob("**/*.csv"), key=lambda p: p.stat().st_mtime, reverse=True))
        if candidates:
            out_in_run = candidates[0]
            session.log(f"Found alternative output: {out_in_run}")
        else:
            session.error("❌ No CSV output found at all inside run_dir.")

    # Print output content
    if out_in_run.exists():
        session.log("\n=== OUTPUT FILE CONTENT ===")
        session.run("bash", "-c", f"echo '--- {out_in_run} ---'; cat {shlex.quote(str(out_in_run))}", external=True)

        out_final.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(out_in_run, out_final)
        session.log(f"\n✅ Copied output to: {out_final}")
    else:
        session.error(f"❌ Output file still missing: {out_in_run}")

# Utilities needed at top level
import shutil  # after function defs to keep nox import first
