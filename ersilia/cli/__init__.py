# Process-wide settings that only the CLI needs. They are applied here, at the
# CLI entry, rather than when the ersilia package is imported, so that using
# ersilia as a library (ersilia.api) does not change the caller's process:
# - session and model files are shared with model containers, so they are
#   created without permission restrictions;
# - ersilia itself runs on CPU;
# - the shell profile gets the CLI snippet.
import os

os.umask(0)
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

from ..default import bashrc_cli_snippet  # noqa: E402
from ..utils.session import create_session_dir  # noqa: E402

bashrc_cli_snippet(overwrite=False)
from .create_cli import create_ersilia_cli  # noqa: E402
from .echo import echo, spinner  # noqa: E402

cli = create_ersilia_cli()
create_session_dir()
if __name__ == "__main__":
    cli()

__all__ = ["echo", "spinner"]
