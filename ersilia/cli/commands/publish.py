import click

from . import ersilia_cli
from .. import echo
from ...publish.publish import ModelPublisher
from ...publish.lake import LakeStorer
from ... import ModelBase


def publish_cmd():
    """Create publish commmand"""

    def _publish(mf, model_id):
        mf.publish(model_id)

    # Example usage: ersilia publish {STEP} {MODEL}
    @ersilia_cli.command(
        short_help="Publish model", help="Contribute a model to the Ersilia Model Hub"
    )
    @click.argument("step", type=click.STRING)
    @click.argument("model", type=click.STRING)
    def publish(step, model):
        model_id = ModelBase(model).model_id
        mp = ModelPublisher(model_id, config_json=None, credentials_json=None)
        ls = LakeStorer(model_id, config_json=None, credentials_json=None)
        if step == "create":
            mp.create()
        elif step == "rebase":
            mp.rebase()
        elif step == "push":
            mp.push()
        elif step == "store":
            ls.store()
        elif step == "test":
            mp.test()
        else:
            echo(
                "Step {0} is not valid. Please choose one of 'create', 'rebase', 'push', 'store' and 'test'",
                fg="red",
            )
        echo(
            ":thumbs_up: Publishing step {0} for model {1} done successfully!".format(
                step, model_id
            ),
            fg="green",
        )
