from .create_cli import create_ersilia_cli
from .echo import echo
import subprocess as s
import sys

cli = create_ersilia_cli()

# Check if Git LFS is installed
def check_gitlfs_installed():
    try:
        check = s.run(["git-lfs"], capture_output=True)
        
    except ModuleNotFoundError:
        sys.exit('\nGit LFS is not installed. We recommend installing Git LFS to'\
            ' use large size models.\n\nCheck out https://git-lfs.github.com/'\
            ' on how to install. After installation, simply use the command'\
            ' `git lfs install` in your repository.\n')

check_gitlfs_installed()

if __name__ == "__main__":
    cli()
