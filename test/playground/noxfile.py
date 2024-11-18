import nox
import yaml
import shutil
from pathlib import Path
from ersilia.utils.logging import logger

ORIGINAL_DIR = Path.cwd()  
config_path = Path("config.yml")
config = yaml.safe_load(config_path.read_text())
REPO_URL = "https://github.com/ersilia-os/ersilia.git"
REPO_DIR = Path("ersilia")

def update_yaml_values(new_values: dict):
    config = yaml.safe_load(config_path.read_text())
    config.update(new_values)
    config_path.write_text(yaml.dump(config))
    logger.info(f"Updated config.yml with: {new_values}")

@nox.session(venv_backend="conda")
def test_cli(session):
    session.install(
        "pytest", 
        "pytest-asyncio", 
        "pytest-xdist", 
        "psutil",
        "PyYAML",
        "rich"
    )

    if config.get("use_existing_env", False):
        logger.info("Using existing environment, skipping setup and installation.")
        
        session.run("ersilia", "--help", external=True)
        session.run(
            "pytest", 
            "commands.py", 
            "-v", 
            silent=False,
            external=True
        )
        return

    session.conda_install(f"python={config.get('python_version', '3.10')}")

    if REPO_DIR.exists() and config.get("overwrite_ersilia_repo", False):
        logger.info(f"Overwriting existing repository directory: {REPO_DIR}")
        shutil.rmtree(REPO_DIR)

    if not REPO_DIR.exists():
        session.run(
            "git", 
            "clone", 
            REPO_URL, 
            external=True
        )
    else:
        logger.info(f"Using existing repository directory: {REPO_DIR}")
    
    session.chdir(REPO_DIR)
    session.install("-e", ".")
    session.chdir(ORIGINAL_DIR)

    session.run(
        "ersilia", 
        "--help", 
        external=True
    )
    # default run --from_github
    logger.info(f"CLI test for model: {config.get("model_id")} and {config.get("fetch_flags")}")
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )
    # default run --from_dockerhub
    update_yaml_values({
        "fetch_flags": "--from_dockerhub"
    })

    logger.info(f"CLI test for model: {config.get("model_id")} and {config.get("fetch_flags")}")
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )
    # CLI test run for auto fetcher decider
    update_yaml_values({
        "fetch_flags": ""
    })
    
    logger.info(f"CLI test for model: {config.get("model_id")} and auto fetcher decider")
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )

    logger.info(f"Fetching and Serving Multiple Models: Fetching")
    update_yaml_values({
        "runner": "multiple",
        "cli_type": "fetch",
        "fetch_flags": "--from_dockerhub"
    })
    
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )

    logger.info(f"Fetching and Serving Multiple Models: Serving")
    update_yaml_values({
        "runner": "multiple",
        "cli_type": "serve"
    })
    
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )

    logger.info(f"Standard and Conventional Run: Conventional")
    update_yaml_values({
        "runner": "single",
        "cli_type": "all",
        "fetch_flags": "--from_dockerhub",
        "output_file": "files/output_eos9gg2_0.json",
        "output_redirection": "true"
    })
    
    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )