import os
import shutil
import sys
import warnings

warnings.filterwarnings("ignore", message="Using slow pure-python SequenceMatcher")

# ruff: noqa
MISSING_PACKAGES = False
try:
    from fuzzywuzzy import fuzz
    from rich.console import Console
    from rich.table import Table
    from rich.text import Text
    from scipy.stats import spearmanr
except ImportError:
    MISSING_PACKAGES = True
# ruff: enable

from ... import ErsiliaBase
from ...default import (
    EOS_TMP,
)

from .services.setup import SetupService
from .services.inspect import InspectService
from .services.io import IOService
from .services.checks import CheckService
from .services.runner import RunnerService
from ...core.modelbase import ModelBase
from ...cli import echo


class ModelTester(ErsiliaBase):
    """
    Class to handle model testing. Initializes the model tester services and runs the tests.

    Parameters
    ----------
    model : str
        The ID of the model.
    from_dir : str
        The directory for the model.
    from_github : bool
        Flag indicating whether to fetch the repository from GitHub.
    from_dockerhub : bool
        Flag indicating whether to fetch the repository from DockerHub.
    from_s3 : bool
        Flag indicating whether to fetch the repository from S3.
    version : str
        Version of the model.
    shallow : bool
        Flag indicating whether to perform shallow checks.
    deep : bool
        Flag indicating whether to perform deep checks.
    surface : bool
        Flag indicating whether to perform surface checks.
    inspect : bool
        Flag indicating whether to perform inspect checks.
    report_path : str
        Flag to specify the path for output as json.
    clean : bool
        Flag indicating whether to clean temp folder.
    """

    def __init__(
        self,
        model,
        from_dir,
        from_github,
        from_dockerhub,
        from_s3,
        version,
        shallow,
        deep,
        surface,
        inspect,
        report_path,
        clean,
    ):
        ErsiliaBase.__init__(self, config_json=None, credentials_json=None)
        clean and self.clean_temp()
        ModelBase(model_id_or_slug=model).is_valid()
        self.model_id = model
        self.from_dir = from_dir
        self.model_dir = os.path.join(EOS_TMP, self.model_id)
        self.dir = from_dir or self.model_dir
        self.default_source = not any([from_dir, from_dockerhub, from_github, from_s3])
        self.from_github = True if self.default_source else from_github
        self.from_dockerhub = from_dockerhub
        self.from_s3 = from_s3
        self.version = version
        self.shallow = shallow
        self.deep = deep
        self.surface = surface
        self.inspect = inspect
        self.report_path = report_path
        self._check_pedendency()
        self.setup_service = SetupService(
            self.model_id,
            self.dir,
            self.from_github,
            self.from_s3,
            self.logger,
        )
        self.ios = IOService(self.logger, self.model_id, self.dir, self.from_dir)
        self.checks = CheckService(
            self.logger,
            self.model_id,
            self.dir,
            self.from_github,
            self.from_s3,
            self.ios,
        )
        self.inspector = InspectService(dir=self.dir, model=self.model_id, remote=True)
        self.runner = RunnerService(
            self.model_id,
            self.logger,
            self.ios,
            self.checks,
            self.setup_service,
            self.dir,
            self.from_dir,
            self._model_path,
            self.from_github,
            self.from_s3,
            self.from_dockerhub,
            self.version,
            self.shallow,
            self.deep,
            self.report_path,
            self.inspector,
            self.surface,
            self.inspect
        )

    def _check_pedendency(self):
        if MISSING_PACKAGES:
            raise ImportError(
                "Missing packages required for testing. "
                "Please install test extras with 'pip install ersilia[test]'."
            )

    def run(self):
        """
        Run the model tester.
        """
        self.runner.run()

    def clean_temp(self):
        """
        Clean up the temp model directory.
        """
        echo("cleaning the temp directory...", fg="yellow", bold=True)
        if os.path.exists(EOS_TMP):
            shutil.rmtree(EOS_TMP)
        else:
            echo("Temp directory does not exist!", fg="yellow", bold=True)
        sys.exit(
            echo(
                "Test command existed after cleaning temp directory.",
                fg="green",
                bold=True,
            )
        )
