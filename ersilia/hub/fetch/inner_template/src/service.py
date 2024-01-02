from typing import List

from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.types import JsonSerializable
from bentoml.service import BentoServiceArtifact

import pickle
import os
import shutil
import tempfile
import subprocess
import csv

CHECKPOINTS_BASEDIR = "checkpoints"
FRAMEWORK_BASEDIR = "framework"


def load_model(framework_dir, checkpoints_dir):
    mdl = Model()
    mdl.load(framework_dir, checkpoints_dir)
    return mdl


def Float(x):
    try:
        return float(x)
    except:
        return None


def String(x):
    x = str(x)
    if not x:
        return None
    if x == "nan":
        return None
    if x == "null":
        return None
    if x == "False":
        return None
    if x == "None":
        return None
    return x


class Model(object):
    def __init__(self):
        self.DATA_FILE = "_data.csv"
        self.OUTPUT_FILE = "_output.csv"
        self.RUN_FILE = "_run.sh"
        self.LOG_FILE = "run.log"

    def load(self, framework_dir, checkpoints_dir):
        self.framework_dir = framework_dir
        self.checkpoints_dir = checkpoints_dir

    def set_checkpoints_dir(self, dest):
        self.checkpoints_dir = os.path.abspath(dest)

    def set_framework_dir(self, dest):
        self.framework_dir = os.path.abspath(dest)

    def run(self, input_list):
        tmp_folder = tempfile.mkdtemp(prefix="eos-")
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
    def __init__(self, name):
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

    def pack(self, model):
        self._model = model
        return self

    def load(self, path):
        model_file_path = self._model_file_path(path)
        model = pickle.load(open(model_file_path, "rb"))
        model.set_checkpoints_dir(
            os.path.join(os.path.dirname(model_file_path), "checkpoints")
        )
        model.set_framework_dir(
            os.path.join(os.path.dirname(model_file_path), "framework")
        )
        return self.pack(model)

    def get(self):
        return self._model

    def save(self, dst):
        self._copy_checkpoints(dst)
        self._copy_framework(dst)
        pickle.dump(self._model, open(self._model_file_path(dst), "wb"))


@artifacts([Artifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def run(self, input: List[JsonSerializable]):
        input = input[0]
        input_list = [inp["input"] for inp in input]
        output = self.artifacts.model.run(input_list)
        return [output]
