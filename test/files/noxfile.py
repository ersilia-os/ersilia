import nox
import subprocess
import os
import re
from pathlib import Path

root = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(root, "..", "..", "README.md")


def extract_shell_commands(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    shell_code_blocks = re.findall(
        r"```(?:bash|sh|shell)\n(.*?)```", content, re.DOTALL
    )

    all_commands = []

    for block in shell_code_blocks:
        commands = [
            cmd for cmd in block.strip().splitlines() if cmd and not cmd.startswith("#")
        ]
        all_commands.extend(commands)

    return " && ".join(all_commands)


@nox.session(venv_backend="conda")
def test_readme_simple(session):
    readme_file_path = readme_path
    commands = extract_shell_commands(readme_file_path)

    result = subprocess.run(
        commands,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Command failed: {commands}\nError: {result.stderr}"


@nox.session(venv_backend="conda")
def test_ersilia_book(session):
    # session.conda_install("python=3.10")
    # session.install("pytest")
    # session.run(
    #     "git", "clone", "https://github.com/ersilia-os/ersilia.git", external=True
    # )
    # session.chdir("ersilia")
    # session.install("-e", ".")
    # session.install("isaura==0.1")
    session.run("pytest", "test_gitbook.py")
