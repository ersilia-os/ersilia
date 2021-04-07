import os
import sys
import click

from ..default import EOS
from ..auth.auth import Auth

from .cmd import Command
from .commands import ersilia_cli


def create_ersilia_cli(pip_installed_bundle_path=None):
    is_contributor = Auth().is_contributor()
    is_contributor = False

    cmd = Command()

    cmd.auth()
    cmd.fetch()
    cmd.delete()
    cmd.serve()
    cmd.close()
    cmd.api()
    cmd.catalog()
    cmd.card()

    # TODO: check if conda and docker are installed.
    if is_contributor:
        cmd.conda()
        cmd.dockerize()

    # Functions only for contributors
    if is_contributor:
        cmd.store()
        cmd.deploy()
        cmd.setup()

    return ersilia_cli
