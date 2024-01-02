from ... import ErsiliaBase
import os
import tempfile
from ...utils.terminal import run_command


class ErsiliaError(Exception):
    """Base class for managing errors in Ersilia"""

    def __init__(
        self, message="Ersilia has experienced an error", hints="", config_json=None
    ):
        text = "Ersilia exception class:\n"
        text += "{}\n\n".format(type(self).__name__)
        text += "Detailed error:\n"
        text += message
        text += "\n\n"
        if hints:
            text += "Hints:\n"
            text += hints
            text += "\n"
        eb = ErsiliaBase(config_json=config_json, credentials_json=None)
        eb.logger.error(text)
        Exception.__init__(self, text)


class MissingDependencyError(ErsiliaError):
    def __init__(self, dependency):
        self.dependency = dependency
        self.message = "Missing dependency {0}".format(self.dependency)
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class NullModelIdentifierError(ErsiliaError):
    def __init__(self, model):
        self.model = model
        self.message = "Model identifier {0} is null".format(self.model)
        self.hints = "This type of error typically occurs when a model has not been served. Please run 'ersilia serve MODEL_ID' if you have a model identifier in mind"
        ErsiliaError.__init__(self, self.message, self.hints)


class InvalidModelIdentifierError(ErsiliaError):
    def __init__(self, model):
        self.model = model
        self.message = "Could not identify model identifier or slug: {0}:".format(
            self.model
        )
        self.hints = "Please check that {0} exists in the Ersilia Model Hub:\n - https://ersilia.io/model-hub (for approved models)\n - https://airtable.com/shrUcrUnd7jB9ChZV (for approved and in preparation models)".format(
            self.model
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelNotAvailableLocallyError(ErsiliaError):
    def __init__(self, model):
        self.model = model
        self.message = (
            "Model {0} is not available locally, so it cannot be served".format(
                self.model
            )
        )
        self.hints = "Fetch the model using the CLI. Simply run:\n"
        self.hints += "$ ersilia fetch {0}".format(self.model)
        ErsiliaError.__init__(self, self.message, self.hints)


class EmptyOutputError(ErsiliaError):
    def __init__(self, model_id, api_name):
        self.model_id = model_id
        self.api_name = api_name
        self.message = "Model API {0}:{1} did not produce an output".format(
            self.model_id, self.api_name
        )
        log = self.run_from_terminal()
        self.message += log
        self.hints = "- Visit the fetch troubleshooting site"
        ErsiliaError.__init__(self, self.message, self.hints)

    def run_from_terminal(self):
        eb = ErsiliaBase()
        bundle_dir = eb._get_bundle_location(model_id=self.model_id)
        framework_dir = os.path.join(
            bundle_dir, self.model_id, "artifacts", "framework"
        )
        bash_executables = ["run.sh", "run_predict.sh", "run_calculate.sh"]
        for exec_file in os.listdir(framework_dir):
            if exec_file in bash_executables:
                break
        exec_file = os.path.join(framework_dir, exec_file)
        input_file = os.path.join(framework_dir, "example_input.csv")
        output_file = os.path.join(framework_dir, "example_output.csv")
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        log_file = os.path.join(tmp_folder, "terminal.log")
        run_command("ersilia example {0} -n 3 -f {1}".format(self.model_id, input_file))
        cmd = "bash {0} {1} {2} {3} 2>&1 | tee -a {4}".format(
            exec_file, framework_dir, input_file, output_file, log_file
        )
        run_command(cmd)
        with open(log_file, "r") as f:
            log = f.read()
        return log
