from ..auth.auth import Auth
from .cmd import Command
from .commands import ersilia_cli


def create_ersilia_cli():
    is_contributor = Auth().is_contributor()

    cmd = Command()

    cmd.auth()
    cmd.catalog()
    cmd.clear()
    cmd.close()
    cmd.delete()
    cmd.example()
    cmd.fetch()
    cmd.info()
    cmd.test()
    cmd.inspect()

    # TODO: publishing functionalities
    if is_contributor:
        cmd.publish()

    cmd.sample()
    cmd.serve()
    cmd.run()

    # TODO: functions only for contributors
    # Functions only for contributors
    if is_contributor:
        cmd.setup()

    return ersilia_cli
