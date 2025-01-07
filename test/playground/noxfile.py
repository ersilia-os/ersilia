import shutil
from pathlib import Path

import nox
import yaml

from ersilia.utils.logging import logger

ORIGINAL_DIR = Path.cwd()
config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())
REPO_URL = "https://github.com/ersilia-os/ersilia.git"
REPO_DIR = Path("ersilia")


def update_yaml_values(new_values: dict):
    existing_config = yaml.safe_load(config_path.read_text())
    existing_config.update(new_values)
    config_path.write_text(yaml.dump(existing_config))


def get_python_version():
    return config.get("python_version", "3.10.10")


def install_dependencies(session):
    session.install(
        "pytest",
        "pytest-asyncio",
        "pytest-xdist",
        "psutil",
        "PyYAML",
        "rich",
    )


def setup_ersilia(session):
    if REPO_DIR.exists() and config.get("overwrite_ersilia_repo", False):
        logger.info(f"Overwriting existing repository directory: {REPO_DIR}")
        shutil.rmtree(REPO_DIR)

    if not REPO_DIR.exists():
        session.run(
            "git",
            "clone",
            REPO_URL,
            external=True,
        )
    else:
        logger.info(f"Using existing repository directory: {REPO_DIR}")

    session.chdir(REPO_DIR)
    session.install("-e", ".")
    session.chdir(ORIGINAL_DIR)


def parse_and_update_config(session):
    passed_args = {
        arg.split("=")[0]: arg.split("=")[1] for arg in session.posargs if "=" in arg
    }

    for key, value in passed_args.items():
        if value.lower() in ("true", "false"):
            passed_args[key] = value.lower() == "true"
        elif "," in value:
            passed_args[key] = value.split(",")

    new_config = {key: passed_args.get(key, config.get(key)) for key in config.keys()}

    logger.info("Parsed arguments and updated configuration:")
    for key in passed_args.keys():
        logger.info(f"  {key}: {passed_args[key]}")

    update_yaml_values(new_config)

    return new_config


@nox.session(venv_backend="conda", python=get_python_version(), reuse_venv=True)
def setup(session):
    install_dependencies(session)
    setup_ersilia(session)


@nox.session(venv_backend="conda", python=get_python_version())
def test_from_github(session):
    install_dependencies(session)
    setup_ersilia(session)
    update_yaml_values({"fetch_flags": "--from_github"})
    logger.info(
        f'CLI test for model: {config.get("model_id")} and {config.get("fetch_flags")}'
    )
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)


@nox.session(venv_backend="conda", python=get_python_version(), reuse_venv=True)
def test_from_dockerhub(session):
    install_dependencies(session)
    setup_ersilia(session)
    update_yaml_values({"fetch_flags": "--from_dockerhub"})
    logger.info(
        f'CLI test for model: {config.get("model_id")} and {config.get("fetch_flags")}'
    )
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)


@nox.session(venv_backend="conda", python=get_python_version())
def test_auto_fetcher_decider(session):
    install_dependencies(session)
    setup_ersilia(session)
    update_yaml_values({"fetch_flags": ""})
    logger.info(
        f'CLI test for model: {config.get("model_id")} and auto fetcher decider'
    )
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)


@nox.session(venv_backend="conda", python=get_python_version())
def test_fetch_multiple_models(session):
    update_yaml_values(
        {
            "runner": "multiple",
            "cli_type": "fetch",
            "fetch_flags": "--from_dockerhub",
        }
    )
    logger.info("Fetching and Serving Multiple Models: Fetching")
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)


@nox.session(venv_backend="conda", python=get_python_version())
def test_serve_multiple_models(session):
    update_yaml_values(
        {"runner": "multiple", "cli_type": "serve", "delete_model": False}
    )
    logger.info("Fetching and Serving Multiple Models: Serving")
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)


@nox.session(venv_backend="conda", python=get_python_version())
def test_conventional_run(session):
    update_yaml_values(
        {
            "runner": "single",
            "cli_type": "all",
            "fetch_flags": "--from_dockerhub",
            "output_files": [
                "files/output_eos9gg2_0.json",
                "files/output_eos9gg2_1.h5",
            ],
            "delete_model": True,
        }
    )
    logger.info("Standard and Conventional Run: Conventional")
    parse_and_update_config(session)
    session.run("pytest", "commands.py", "-v", silent=False)
