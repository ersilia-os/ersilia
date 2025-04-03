import click

from ... import ModelBase
from ...publish.publish import ModelPublisher
from .. import echo
from . import ersilia_cli


def publish_cmd():
    """
    Publishes a specified model.

    This command allows users to publish a specified model to the model hub.

    Returns
    -------
    function
        The publish command function to be used by the CLI.

    Examples
    --------
    .. code-block:: console

        Publish a model by its ID:
        $ ersilia publish create <model_id>

        Rebase a model:
        $ ersilia publish rebase <model_id>
    """

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
        if step == "create":
            mp.create()
        elif step == "rebase":
            mp.rebase()
        elif step == "push":
            mp.push()
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
