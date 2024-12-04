# import subprocess
# import pytest
# import re


# def extract_shell_commands(file_path):
#     with open(file_path, "r") as f:
#         content = f.read()
#     shell_code_blocks = re.findall(
#         r"```(?:bash|sh|shell)\n(.*?)```", content, re.DOTALL
#     )
#     shell_commands = []
#     for block in shell_code_blocks:
#         commands = [
#             cmd for cmd in block.strip().splitlines() if cmd and not cmd.startswith("#")
#         ]
#         shell_commands.extend(commands)
#     return shell_commands


# @pytest.fixture
# def readme_file_path():
#     return "../../README.md"


# def test_shell_commands_in_readme(readme_file_path):
#     commands = extract_shell_commands(readme_file_path)

#     for cmd in commands:
#         if "conda" in cmd:
#             print(f"Skipping Conda command: {cmd}")
#             continue

#         result = subprocess.run(cmd.strip(), shell=True, capture_output=True, text=True)

#         assert result.returncode == 0, f"Command failed: {cmd}\nError: {result.stderr}"
