from ..auth.auth import Auth
from .cmd import Command
from .commands import ersilia_cli


def create_ersilia_cli(pip_installed_bundle_path=None):
    is_contributor = Auth().is_contributor()

    cmd = Command()

    cmd.auth()
    cmd.fetch()
    cmd.delete()
    cmd.serve()
    cmd.close()
    cmd.api()
    cmd.example()
    cmd.catalog()
    cmd.card()

    if is_contributor:
        cmd.publish()

    # TODO: check if conda and docker are installed.
    if is_contributor:
        cmd.conda()
        cmd.dockerize()

    # TODO: functions only for contributors
    # Functions only for contributors
    if is_contributor:
        cmd.setup()

    return ersilia_cli
