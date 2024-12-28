from ..utils.session import create_session_dir
from .create_cli import create_ersilia_cli
from .echo import echo

cli = create_ersilia_cli()
create_session_dir()
if __name__ == "__main__":
    cli()

__all__ = ["echo"]
