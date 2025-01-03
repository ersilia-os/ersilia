import csv
import os
import pickle
import shutil
import subprocess
from typing import Any, Dict, List

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service import BentoServiceArtifact
from bentoml.types import JsonSerializable

from .....utils.logging import make_temp_dir

CHECKPOINTS_BASEDIR = "checkpoints"
FRAMEWORK_BASEDIR = "framework"


def load_model(framework_dir: str, checkpoints_dir: str) -> "Model":
    """
    Load the model with the given framework and checkpoints directories.

    Parameters
    ----------
    framework_dir : str
        Path to the framework directory.
    checkpoints_dir : str
        Path to the checkpoints directory.

    Returns
    -------
    Model
        Loaded model instance.
    """
    mdl = Model()
    mdl.load(framework_dir, checkpoints_dir)
    return mdl


def Float(x: Any) -> float:
    """
    Convert the input to a float.

    Parameters
    ----------
    x : Any
        Input value to be converted.

    Returns
    -------
    float
        Converted float value or None if conversion fails.
    """
    try:
        return float(x)
    except:
        return None


def String(x: Any) -> str:
    """
    Convert the input to a string.

    Parameters
    ----------
    x : Any
        Input value to be converted.

    Returns
    -------
    str
        Converted string value or None if conversion fails.
    """
    x = str(x)
    if not x or x in ["nan", "null", "False", "None"]:
        return None
    return x


class Model(object):
    """
    A class used to represent the Model.

    Attributes
    ----------
    DATA_FILE : str
        Filename for the data file.
    OUTPUT_FILE : str
        Filename for the output file.
    RUN_FILE : str
        Filename for the run script.
    LOG_FILE : str
        Filename for the log file.

    Methods
    -------
    load(framework_dir, checkpoints_dir)
        Loads the model with the given framework and checkpoints directories.
    set_checkpoints_dir(dest)
        Sets the checkpoints directory.
    set_framework_dir(dest)
        Sets the framework directory.
    run(input_list)
        Runs the model with the given input list and returns the result.
    """

    def __init__(self):
        self.DATA_FILE = "_data.csv"
        self.OUTPUT_FILE = "_output.csv"
        self.RUN_FILE = "_run.sh"
        self.LOG_FILE = "run.log"

    def load(self, framework_dir: str, checkpoints_dir: str):
        """
        Load the model with the given framework and checkpoints directories.

        Parameters
        ----------
        framework_dir : str
            Path to the framework directory.
        checkpoints_dir : str
            Path to the checkpoints directory.
        """
        self.framework_dir = framework_dir
        self.checkpoints_dir = checkpoints_dir

    def set_checkpoints_dir(self, dest: str):
        """
        Set the checkpoints directory.

        Parameters
        ----------
        dest : str
            Path to the checkpoints directory.
        """
        self.checkpoints_dir = os.path.abspath(dest)

    def set_framework_dir(self, dest: str):
        """
        Set the framework directory.

        Parameters
        ----------
        dest : str
            Path to the framework directory.
        """
        self.framework_dir = os.path.abspath(dest)

    def run(self, input_list: List[str]) -> Dict[str, Any]:
        """
        Run the model with the given input list and return the result.

        Parameters
        ----------
        input_list : List[str]
            List of input strings to be processed by the model.

        Returns
        -------
        dict
            Dictionary containing the result and metadata.
        """
        tmp_folder = make_temp_dir(prefix="eos-")
        data_file = os.path.join(tmp_folder, self.DATA_FILE)
        output_file = os.path.join(tmp_folder, self.OUTPUT_FILE)
        log_file = os.path.join(tmp_folder, self.LOG_FILE)
        with open(data_file, "w") as f:
            f.write("input" + os.linesep)
            for inp in input_list:
                f.write(inp + os.linesep)
        run_file = os.path.join(tmp_folder, self.RUN_FILE)
        with open(run_file, "w") as f:
            lines = [
                "bash {0}/run.sh {0} {1} {2}".format(
                    self.framework_dir, data_file, output_file
                )
            ]
            f.write(os.linesep.join(lines))
        cmd = "bash {0}".format(run_file)
        with open(log_file, "w") as fp:
            subprocess.Popen(
                cmd, stdout=fp, stderr=fp, shell=True, env=os.environ
            ).wait()
        with open(output_file, "r") as f:
            reader = csv.reader(f)
            h = next(reader)
            R = []
            for r in reader:
                R += [
                    {"outcome": [Float(x) for x in r]}
                ]  # <-- EDIT: Modify according to type of output (Float, String...)
        meta = {"outcome": h}
        result = {"result": R, "meta": meta}
        shutil.rmtree(tmp_folder)
        return result


class Artifact(BentoServiceArtifact):
    """
    A class used to represent the Artifact.

    Attributes
    ----------
    name : str
        Name of the artifact.
    _model : Model
        Model instance.
    _extension : str
        File extension for the artifact.

    Methods
    -------
    pack(model)
        Packs the model into the artifact.
    load(path)
        Loads the model from the given path.
    get()
        Returns the model instance.
    save(dst)
        Saves the model to the given destination.
    """

    def __init__(self, name: str):
        BentoServiceArtifact.__init__(name)
        self._model = None
        self._extension = ".pkl"

    def _copy_checkpoints(self, base_path):
        src_folder = self._model.checkpoints_dir
        dst_folder = os.path.join(base_path, "checkpoints")
        if os.path.exists(dst_folder):
            os.rmdir(dst_folder)
        shutil.copytree(src_folder, dst_folder)

    def _copy_framework(self, base_path):
        src_folder = self._model.framework_dir
        dst_folder = os.path.join(base_path, "framework")
        if os.path.exists(dst_folder):
            os.rmdir(dst_folder)
        shutil.copytree(src_folder, dst_folder)

    def _model_file_path(self, base_path):
        return os.path.join(base_path, self.name + self._extension)

    def pack(self, model: Model) -> "Artifact":
        """
        Pack the model into the artifact.

        Parameters
        ----------
        model : Model
            Model instance to be packed.

        Returns
        -------
        Artifact
            The artifact instance with the packed model.
        """
        self._model = model
        return self

    def load(self, path: str) -> "Artifact":
        """
        Load the model from the given path.

        Parameters
        ----------
        path : str
            Path to load the model from.

        Returns
        -------
        Artifact
            The artifact instance with the loaded model.
        """
        model_file_path = self._model_file_path(path)
        model = pickle.load(open(model_file_path, "rb"))
        model.set_checkpoints_dir(
            os.path.join(os.path.dirname(model_file_path), "checkpoints")
        )
        model.set_framework_dir(
            os.path.join(os.path.dirname(model_file_path), "framework")
        )
        return self.pack(model)

    def get(self) -> Model:
        """
        Get the model instance.

        Returns
        -------
        Model
            The model instance.
        """
        return self._model

    def save(self, dst: str):
        """
        Save the model to the given destination.

        Parameters
        ----------
        dst : str
            Destination path to save the model.
        """
        self._copy_checkpoints(dst)
        self._copy_framework(dst)
        pickle.dump(self._model, open(self._model_file_path(dst), "wb"))


@artifacts([Artifact("model")])
class Service(BentoService):
    """
    A class used to represent the Service.

    Methods
    -------
    run(input)
        Runs the service with the given input and returns the output.
    """

    @api(input=JsonInput(), batch=True)
    def run(self, input: List[JsonSerializable]) -> List[Dict[str, Any]]:
        """
        Run the service with the given input and return the output.

        Parameters
        ----------
        input : List[JsonSerializable]
            List of JSON serializable input data.

        Returns
        -------
        List[Dict[str, Any]]
            List containing the output data.
        """
        input = input[0]
        input_list = [inp["input"] for inp in input]
        output = self.artifacts.model.run(input_list)
        return [output]
