from ersilia import __version__
import sys
import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from bentoml.server import start_dev_server
from bentoml.cli.click_utils import conditional_argument
from bentoml.cli.bento_service import resolve_bundle_path
from bentoml.saved_bundle import load_bento_service_api


def create_ersilia_service_cli(pip_installed_bundle_path=None):
    from ersilia.auth.auth import Auth

    is_contributor = Auth().is_contributor()

    @click.group(cls=BentoMLCommandGroup)
    @click.version_option(version=__version__)
    def ersilia_cli():
        """
        Ersilia CLI tool
        """

    # Example usage: ersilia auth login
    @ersilia_cli.command(
        short_help="Log in to ersilia to enter contributor mode",
        help="Log in to ersilia to enter contributor mode. "
             "GitHub credentials are used."
    )
    def auth():
        pass

    # Example usage: ersilia conda activate eos0aaa
    @ersilia_cli.command(
        short_help="Use conda from ersilia.",
        help="Use conda from ersilia in order to access model environments."
    )
    def conda():
        pass

    # Example usage: ersilia fetch {MODEL_ID}
    @ersilia_cli.command(
        short_help="Fetch model from EOS repository",
        help="Fetch model from EOS repository. Model files are downloaded from GitHub and model data are "
             "downloaded from a file storage system such as the Open Science Framework. Model is downloaded to "
             "an EOS folder, then packed to a BentoML bundle",
    )
    @click.argument("model_id", type=click.STRING)
    def fetch(model_id):
        from ersilia.hub.fetch import ModelFetcher
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
        from ersilia.hub.delete import ModelFullDeleter
        click.echo(click.style("Deleting model %s" % model_id, fg="yellow"))
        ModelFullDeleter().delete(model_id)

    # Example usage: ersilia serve {MODEL_ID}
    @ersilia_cli.command(
        help="Serve model using BentoML functionality",
        short_help="Serve model",
    )
    @conditional_argument(pip_installed_bundle_path is None, "model_id", type=click.STRING)
    @click.option(
        "--port",
        type=click.INT,
        default=None,
        help="The port to listen on for the REST api server, "
             "default is to look for a free port",
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
        if port is None:
            from ersilia.utils.ports import find_free_port
            port = find_free_port()
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
    @click.option(
        '--hub',
        is_flag=True,
        default=True,
        help="Show catalog of models available in the Ersilia Hub"
    )
    def catalog(local=False, hub=True):
        from ersilia.hub.catalog import ModelCatalog
        mc = ModelCatalog(as_dataframe=False)
        if not local:
            click.echo(mc.hub())
        else:
            click.echo(mc.local())

    # Example usage: ersilia card {MODEL_ID}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia hub."
    )
    @click.argument("model_id", type=click.STRING)
    def card(model_id):
        from ersilia.hub.card import ModelCard
        mc = ModelCard()
        click.echo(mc.get(model_id, as_json=True))

    # Example usage: ersilia app {MODEL_ID}
    @ersilia_cli.command(
        short_help="Run model app",
        help="Run app designed to run the model in the browser. Apps typically do not allow for large-scale inputs. "
             "Please note that, at the moment, not all models have an app available.",
    )
    @click.argument("model_id", type=click.STRING)
    def app(model_id):
        from ersilia.app.app import StreamlitApp
        sa = StreamlitApp()
        status = sa.run(model_id)
        if not status:
            click.echo(click.style("App could not be run or it is not available for model %s" % model_id, fg="red"))
            click.echo(click.style("Check that an app.py script exists in the model repository", fg="red"))

    # Example usage: ersilia setup
    @ersilia_cli.command(
        short_help="Setup ersilia",
        help="Setup ersilia, including building a model-server image, a base environment (eos), rdkit, etc."
    )
    @click.option(
        '--base',
        is_flag=True,
        default=False,
        help="Install only bare-minimum dependencies."
    )
    @click.option(
        '--full',
        is_flag=True,
        default=True,
        help="Install all the necessary dependencies."
    )
    def setup(base=False, full=True):
        from ersilia.utils.installers import base_installer, full_installer
        if base:
            base_installer()
        elif full:
            full_installer()
        else:
            pass

    # Functions only for contributors
    if is_contributor:
        # Example usage: ersilia store {MODEL_ID}
        @ersilia_cli.command(
            short_help="Store a model",
            help="Store a model in GitHub and the chosen file storage. "
                 "This option is only for developers and requires credentials."
        )
        @click.argument("model_id", type=click.STRING)
        @click.option("--path", required=True, type=click.Path())
        def store(model_id, path):
            from ersilia.contrib.store import ModelStorager
            ms = ModelStorager()
            ms.store(path, model_id)

        # Example usage: ersilia deploy {MODEL_ID}
        @ersilia_cli.command(
            short_help="Deploy model to the cloud",
            help="Deploy model in a cloud service. "
                 "This option is only for developers and requires credentials."
        )
        @click.argument("model_id", type=click.STRING)
        @click.option("--cloud", default="heroku", type=click.STRING)
        def deploy(model_id, cloud):
            from ersilia.contrib.deploy import Deployer
            dp = Deployer(cloud=cloud)
            if dp.dep is None:
                click.echo(click.style("Please enter a valid cloud option", fg="red"))
                click.echo(click.style("Only 'heroku' and 'local' are available for the moment...", fg="yellow"))
                return
            dp.deploy(model_id)

    return ersilia_cli
