from ersilia import __version__
import sys
import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from bentoml.server import start_dev_server
from bentoml.cli.click_utils import conditional_argument
from bentoml.cli.bento_service import resolve_bundle_path
from bentoml.saved_bundle import load_bento_service_api
from ersilia.hub.fetch import ModelFetcher
from ersilia.hub.catalog import ModelCatalog
from ersilia.hub.delete import ModelEosDeleter, ModelTmpDeleter, ModelBentoDeleter, ModelPipDeleter, ModelBundleDeleter
from ersilia.hub.card import ModelCard
from ersilia.app.app import StreamlitApp
from ersilia.core.base import ErsiliaBase
from ersilia.contrib.store import ModelStorager
from ersilia.contrib.deploy import Deployer


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
        click.echo(click.style("Fetching model %s" % model_id, fg="yellow"))
        mf = ModelFetcher()
        mf.fetch(model_id)
        click.echo(click.style("Model fetched successfully!", fg="green"))

    # Example usage: ersilia delete {MODEL_ID}
    @ersilia_cli.command(
        short_help="Delete model from local computer",
        help="Delete model from local computer. The BentoML bundle is deleted, as well as the files stored in "
             "the EOS directory and the Pip-installed package",
    )
    @click.argument("model_id", type=click.STRING)
    def delete(model_id):
        click.echo(click.style("Deleting BentoML files", fg="yellow"))
        ModelBentoDeleter().delete(model_id)
        click.echo(click.style("Deleting EOS files", fg="yellow"))
        ModelEosDeleter().delete(model_id)
        click.echo(click.style("Deleting bundles", fg="yellow"))
        ModelBundleDeleter().delete(model_id)
        click.echo(click.style("Deleting temporary files (if any)", fg="yellow"))
        ModelTmpDeleter().delete(model_id)
        click.echo(click.style("Deleting local Pip package", fg="yellow"))
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
        short_help="Run a prediction",
        help="Run prediction using the BentoML prediction API of the model. See BentoML for more details.",
        context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
    )
    @conditional_argument(pip_installed_bundle_path is None, "model_id", type=click.STRING)
    @click.argument('run_args', nargs=-1, type=click.UNPROCESSED)
    def predict(run_args, model_id=None):
        saved_bundle_path = resolve_bundle_path(model_id + ":latest", pip_installed_bundle_path)
        api = load_bento_service_api(saved_bundle_path, "predict")
        exit_code = api.handle_cli(run_args)
        sys.exit(exit_code)

    # Example usage: ersilia catalog
    @ersilia_cli.command(
        help="List a catalog of models",
    )
    @click.option(
        '--local',
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer"
    )
    def catalog(local=False):
        mc = ModelCatalog(as_dataframe=False)
        if not local:
            click.echo(mc.hub())
        else:
            click.echo(mc.local())

    # Example usage: ersilia app {MODEL_ID}
    @ersilia_cli.command(
        short_help="Run model app",
        help="Run app designed to run the model in the browser. Apps typically do not allow for large-scale inputs. "
             "Please note that, at the moment, not all models have an app available.",
    )
    @click.argument("model_id", type=click.STRING)
    def app(model_id):
        sa = StreamlitApp()
        status = sa.run(model_id)
        if not status:
            click.echo(click.style("App could not be run or it is not available for model %s" % model_id, fg="red"))
            click.echo(click.style("Check that an app.py script exists in the model repository", fg="red"))

    def _is_dev_ready(model_id):
        eb = ErsiliaBase()
        if not eb._has_credentials():
            click.echo(click.style("The store action is reserved to developers... Credentials are needed", fg="red"))
            click.echo(click.style("Please write to us if you want to know more: hello@ersilia.io", fg="blue"))
            return False
        if not eb._is_ready(model_id):
            click.echo(click.style("Model %s does not seem to be ready" % model_id, fg="red"))
            click.echo(click.style("Try running ersilia fetch %s instead" % model_id, fg="blue"))
            return False
        return True

    # Example usage: ersilia store {MODEL_ID}
    @ersilia_cli.command(
        short_help="Store a model [only for developers]",
        help="Store a model in GitHub and the chosen file storage. "
             "This option is only for developers and requires credentials."
    )
    @click.argument("model_id", type=click.STRING)
    @click.option("--path", required=True, type=click.Path())
    def store(model_id, path):
        if not _is_dev_ready(model_id):
            return
        ms = ModelStorager()
        ms.store(path, model_id)

    # Example usage: ersilia deploy {MODEL_ID}
    @ersilia_cli.command(
        short_help="Deploy model to the cloud [only for developers]",
        help="Deploy model in a cloud service. "
             "This option is only for developers and requires credentials."
    )
    @click.argument("model_id", type=click.STRING)
    @click.option("--cloud", default="heroku", type=click.STRING)
    def deploy(model_id, cloud):
        if not _is_dev_ready(model_id):
            return
        dp = Deployer(cloud=cloud)
        if dp.dep is None:
            click.echo(click.style("Please enter a valid cloud option", fg="red"))
            click.echo(click.style("Only 'heroku' and 'local' are available for the moment...", fg="yellow"))
            return
        dp.deploy(model_id)

    # Example usage: ersilia card {MODEL_ID}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia hub."
    )
    @click.argument("model_id", type=click.STRING)
    def card(model_id):
        mc = ModelCard()
        click.echo(mc.get(model_id))

    return ersilia_cli
