from .create_cli import create_ersilia_cli
from .echo import echo
from ..utils.session import create_session_dir

cli = create_ersilia_cli()
create_session_dir()
if __name__ == "__main__":
    cli()
