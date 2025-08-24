import subprocess
from pathlib import Path

import yaml


class NoxSession:
    def __init__(self, name):
        self.name = name

    def execute(self, noxfile):
        try:
            subprocess.run(
                ["nox", "-f", noxfile, "-s", self.name],
                check=True,
            )
            print(f"Session '{self.name}' executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing session '{self.name}': {e}")


class NoxRunner:
    def __init__(self, config_path="config.yml", noxfile="noxfile.py"):
        self.original_dir = Path.cwd()
        self.config_path = Path(config_path)
        self.noxfile = noxfile
        self.config = yaml.safe_load(self.config_path.read_text())
        self.nox_command = "nox"
        self.queue = []

    def update_yaml_values(self, new_values: dict):
        existing_config = yaml.safe_load(self.config_path.read_text())
        existing_config.update(new_values)
        self.config_path.write_text(yaml.dump(existing_config))

    def get_python_version(self):
        return self.config.get("python_version", "3.10.10")

    def add_session(self, session_name):
        self.queue.append(NoxSession(session_name))

    def execute_all(self):
        for session in self.queue:
            session.execute(self.noxfile)
        self.queue.clear()

    def clear_queue(self):
        self.queue.clear()

    def setup(self):
        self.add_session("setup")

    def test_from_github(self):
        self.add_session("test_from_github")

    def test_from_dockerhub(self):
        self.add_session("test_from_dockerhub")

    def test_auto_fetcher_decider(self):
        self.add_session("test_auto_fetcher_decider")

    def test_fetch_multiple_models(self):
        self.add_session("test_fetch_multiple_models")

    def test_serve_multiple_models(self):
        self.add_session("test_serve_multiple_models")
