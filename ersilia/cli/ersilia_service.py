from ersilia import __version__
import sys
import click
import tabulate
from bentoml.cli.click_utils import BentoMLCommandGroup
from bentoml.server import start_dev_server
from bentoml.cli.click_utils import conditional_argument
from bentoml.cli.bento_service import resolve_bundle_path
from bentoml.saved_bundle import load_bento_service_api
from ersilia.hub.fetch import ModelFetcher
from ersilia.hub.list import ModelList
from ersilia.hub.delete import ModelEosDeleter, ModelBentoDeleter, ModelPipDeleter


def create_ersilia_service_cli(pip_installed_bundle_path=None):
    @click.group(cls=BentoMLCommandGroup)
    @click.version_option(version=__version__)
    def ersilia_cli():
        """
        Ersilia CLI tool
        """

    # Example usage: ersilia fetch {MODEL_ID}
    @ersilia_cli.command(
        short_help="Fetch model from EOS repository",
        help="Fetch model from EOS repository. Model files are downloaded from GitHub and model data are "
             "downloaded from a file storage system such as the Open Science Framework. Model is downloaded to "
             "an EOS folder, then packed to a BentoML bundle",
    )
    @click.argument("model_id", type=click.STRING)
    def fetch(model_id):
        mf = ModelFetcher()
        mf.fetch(model_id)

    # Example usage: ersilia delete {MODEL_ID}
    @ersilia_cli.command(
        short_help="Delete model from local computer",
        help="Delete model from local computer. The BentoML bundle is deleted, as well as the files stored in "
             "the EOS directory and the Pip-installed package",
    )
    @click.argument("model_id", type=click.STRING)
    def delete(model_id):
        ModelBentoDeleter().delete(model_id)
        ModelEosDeleter().delete(model_id)
        ModelPipDeleter().delete(model_id)

    # Example usage: ersilia serve {MODEL_ID}
    @ersilia_cli.command(
        help="Serve model using BentoML functionality",
        short_help="Serve model",
    )
    @conditional_argument(pip_installed_bundle_path is None, "model_id", type=click.STRING)
    @click.option(
        "--port",
        type=click.INT,
        default=5000,
        help="The port to listen on for the REST api server, "
             "default is 5000",
        envvar='BENTOML_PORT',
    )
    @click.option(
        '--enable-microbatch/--disable-microbatch',
        default=False,
        help="Run API server with micro-batch enabled",
        envvar='BENTOML_ENABLE_MICROBATCH',
    )
    @click.option(
        '--run-with-ngrok',
        is_flag=True,
        default=False,
        help="Use ngrok to relay traffic on a public endpoint to this "
             "API server on localhost",
        envvar='BENTOML_ENABLE_NGROK',
    )
    def serve(port, model_id=None, enable_microbatch=False, run_with_ngrok=False):
        saved_bundle_path = resolve_bundle_path(model_id + ":latest", pip_installed_bundle_path)
        start_dev_server(saved_bundle_path, port, enable_microbatch, run_with_ngrok)

    # Example usage: ersilia predict --input=...
    @ersilia_cli.command(
        help="Run a prediction",
        short_help="Run prediction",
        context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
    )
    @conditional_argument(pip_installed_bundle_path is None, "model_id", type=click.STRING)
    @click.argument('run_args', nargs=-1, type=click.UNPROCESSED)
    def predict(run_args, model_id=None):
        saved_bundle_path = resolve_bundle_path(model_id + ":latest", pip_installed_bundle_path)
        api = load_bento_service_api(saved_bundle_path, "predict")
        exit_code = api.handle_cli(run_args)
        sys.exit(exit_code)

    # Example usage: ersilia list
    @ersilia_cli.command(
        help="List available models bundled in the BentoML style",
        short_help="List available models"
    )
    def list():
        ml = ModelList()
        df = ml.bentoml()
        R = [["MODEL_ID", "BENTO_SERVICE"]]
        for v in df[["MODEL_ID", "BENTO_SERVICE"]].values:
            R += [[v[0], v[1]]]
        click.echo(tabulate.tabulate(R))

    return ersilia_cli
