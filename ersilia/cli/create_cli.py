from ..auth.auth import Auth
from .cmd import Command
from .commands import ersilia_cli


def create_ersilia_cli():
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
    cmd.clear()

    # TODO: publishing functionalities
    if is_contributor:
        cmd.publish()

    # TODO: functions only for contributors
    # Functions only for contributors
    if is_contributor:
        cmd.setup()

    return ersilia_cli
