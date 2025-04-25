import click

from ...core.session import Session
from ...store.api import InferenceStoreApi
from ...store.utils import OutputSource, echo_intro
from . import ersilia_cli


def dump_cmd():
    """Create dump command"""

    # dumpy usage: ersilia dumpy -n 10 [--output {FILE_NAME}]
    @ersilia_cli.command(
        short_help="Fetch precalculations from cloud and local sources",
        help="This command fetch cached precalculations from a cloud and from local. In the local case precalculation will be fetched from a Redis container.",
    )
    @click.option("--n_samples", "-n", default=-1, type=click.INT)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    def dump(n_samples, output):
        session = Session(config_json=None)
        model_id = session.current_model_id()
        output_source = session.current_output_source()
        echo_intro(
            InferenceStoreApi(click_iface=None, model_id=model_id, output=output).click
        )

        ifst = InferenceStoreApi(
            model_id=model_id,
            output=output,
            output_source=output_source,
            n_samples=n_samples,
        )
        if output_source == OutputSource.CLOUD_ONLY:
            ifst.get_precalculations(None)

    return dump
