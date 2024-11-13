import nox
import yaml
import shutil
from pathlib import Path
from ersilia.utils.logging import logger

ORIGINAL_DIR = Path.cwd()  
config = yaml.safe_load(Path("config.yml").read_text())
REPO_URL = "https://github.com/ersilia-os/ersilia.git"
REPO_DIR = Path("ersilia")

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

    session.run(
        "pytest", 
        "commands.py", 
        "-v", 
        silent=False
    )
