from .create_cli import create_ersilia_cli
from .echo import echo

cli = create_ersilia_cli()

if __name__ == "__main__":
    cli()
