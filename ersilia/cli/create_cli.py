from ersilia import __version__
import os
import sys
import click
from bentoml.cli.click_utils import BentoMLCommandGroup
from bentoml.server import start_dev_server
from bentoml.cli.click_utils import conditional_argument
from bentoml.cli.bento_service import resolve_bundle_path
from bentoml.saved_bundle import load_bento_service_api
from ..default import EOS


def create_ersilia_cli(pip_installed_bundle_path=None):
    from ersilia.auth.auth import Auth

    is_contributor = Auth().is_contributor()
    is_contributor = False

    @click.group(cls=BentoMLCommandGroup)
    @click.version_option(version=__version__)
    def ersilia_cli():
        """
        Ersilia CLI tool
        """

    # Example usage: ersilia auth login
    @ersilia_cli.command(
        short_help="Log in to ersilia to enter contributor mode.",
        help="Log in to ersilia to enter contributor mode. "
             "GitHub credentials are used."
    )
    def auth():
        """
        In the user's system profile there is a redirect to gh auth (GitHub authorisation CLI)
        """
        pass

    if is_contributor:
        # Example usage: ersilia conda activate eos0aaa
        @ersilia_cli.command(
            short_help="Use conda from ersilia.",
            help="Use conda from ersilia in order to access model environments."
        )
        def conda():
            pass

    # Example usage: ersilia fetch {MODEL_ID}
    @ersilia_cli.command(
        short_help="Fetch model from Ersilia Model Hub",
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
        md = ModelFullDeleter()
        if md.needs_delete(model_id):
            click.echo(click.style("Deleting model %s" % model_id, fg="yellow"))
            ModelFullDeleter().delete(model_id)
            click.echo(click.style("Model {0} deleted successfully!".format(model_id), fg="green"))
        else:
            click.echo(click.style("Model {0} is not available locally. No delete is necessary".format(model_id), fg="green"))


    def _tmp_pid_file(model_id):
        return os.path.join(EOS, "tmp", "{0}.pid".format(model_id))

    # Example usage: ersilia serve {MODEL_ID}
    @ersilia_cli.command(
        short_help="Serve model",
        help="Serve model"
    )
    @click.argument("model_id", type=click.STRING)
    def serve(model_id):
        from ersilia.serve.autoservice import AutoService
        import webbrowser
        srv = AutoService(model_id)
        srv.serve()
        click.echo(click.style("Serving model {0}!".format(model_id), fg="green"))
        click.echo()
        click.echo(click.style("    URL: {0}".format(srv.service.url), fg="yellow"))
        click.echo(click.style("    PID: {0}".format(srv.service.pid), fg="yellow"))
        #webbrowser.open(srv.service.url)
        tmp_file = _tmp_pid_file(model_id)
        with open(tmp_file, "a+") as f:
            f.write("{0} {1}".format(srv.service.pid, srv.service.url, os.linesep))

    # Example usage: ersilia close {MODEL_ID}
    @ersilia_cli.command(
        short_help="Close model",
        help="Close model"
    )
    @click.argument("model_id", type=click.STRING)
    def close(model_id):
        from ..utils.terminal import run_command
        tmp_file = _tmp_pid_file(model_id)
        with open(tmp_file, "r") as f:
            for l in f:
                pid = int(l.rstrip().split()[0])
                cmd = "kill {0}".format(pid)
                run_command(cmd, quiet=True)
        os.remove(tmp_file)

    # Example usage: ersilia api {MODEL_ID} {API_NAME} {INPUT}
    @ersilia_cli.command(
        short_help="Run API on a served model",
        help="Run API on a served model"
    )
    @click.argument("model_id", type=click.STRING)
    @click.argument("api_name", type=click.STRING)
    @click.argument("input", type=click.STRING)
    def api(model_id, api_name, input):
        import requests
        tmp_file = _tmp_pid_file(model_id)
        with open(tmp_file, "r") as f:
            for l in f:
                url = l.rstrip().split()[1]
        response = requests.post("{0}/{1}".format(url, api_name), json=input)
        click.echo(response.json())


    # Example usage: ersilia catalog
    @ersilia_cli.command(
        help="List a catalog of models",
    )
    @click.option(
        '--hub',
        is_flag=True,
        default=False,
        help="Show catalog of models available in the Ersilia Model Hub"
    )
    @click.option(
        '--local',
        is_flag=True,
        default=False,
        help="Show catalog of models available in the local computer"
    )
    @click.option(
       '--backlog',
       is_flag=True,
       default=False,
       help="Show models backlog (wish list) in a Google Spreadsheet"
    )
    def catalog(hub=False, local=False, backlog=False):
        from ersilia.hub.catalog import ModelCatalog
        mc = ModelCatalog()
        if hub:
            click.echo(mc.hub())
        elif local:
            click.echo(mc.local())
        elif backlog:
            click.echo(mc.backlog())
        else:
            click.echo(click.style("Specifiy one of --hub --local or --backlog", fg="blue"))

    # Example usage: ersilia card {MODEL_ID}
    @ersilia_cli.command(
        short_help="Get model info card",
        help="Get model info card from Ersilia Model Hub."
    )
    @click.argument("model_id", type=click.STRING)
    def card(model_id):
        from ersilia.hub.card import ModelCard
        mc = ModelCard()
        click.echo(mc.get(model_id, as_json=True))

    if is_contributor:
        # Example usage: ersilia dockerize {MODEL_ID}
        @ersilia_cli.command(
            short_help="Containerize model using docker",
            help="Containerize model in the BentoML style using docker",
        )
        @click.argument("model_id", type=click.STRING)
        def dockerize(model_id):
            from ersilia.hub.fetch import ModelFetcher
            mf = ModelFetcher()
            mf.containerize(model_id=model_id)

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

    return ersilia_cli
