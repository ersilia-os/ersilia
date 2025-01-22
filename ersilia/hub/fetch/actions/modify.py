import os

import yaml

from ersilia.default import PACKMODE_FILE

from ....utils.conda import SimpleConda
from ....utils.logging import make_temp_dir
from ....utils.terminal import run_command
from ...bundle.bundle import (
    BundleDockerfileFile,
    BundleEnvironmentFile,
    BundleRequirementsFile,
)
from .. import DOCKERFILE, ENVIRONMENT_YML
from . import BaseAction


class ModelModifier(BaseAction):
    """
    Modifies model bundles to ensure compatibility and consistency.

    Compatibility: Ensures that the model bundle can run in the target environment without issues.
    Consistency: Ensures that the model bundle follows the same structure and conventions as other bundles.

    Parameters
    ----------
    model_id : str
        Identifier of the model to be modified.
    config_json : dict
        Configuration settings for the modifier.

    Methods
    -------
    modify()
        Modifies the model bundle.
    """

    def __init__(self, model_id: str, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def _bundle_uses_ersilia(self, model_id: str) -> bool:
        src = os.path.join(self._get_bundle_location(model_id), model_id, "src")
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "grep.txt")
        cmd = "grep -R 'ersilia' {0}/* > {1}".format(src, tmp_file)
        run_command(cmd)
        with open(tmp_file, "r") as f:
            grep = f.read()
        if grep:
            return True
        else:
            return False

    def _modify_bundle_environment_yml(self, model_id: str):
        dir = self._get_bundle_location(model_id)
        yml_file = os.path.join(dir, ENVIRONMENT_YML)
        if not os.path.exists(yml_file):
            return
        try:
            with open(yml_file, "r") as f:
                data = yaml.safe_load(f)
                for i, d in enumerate(data["dependencies"][:-1]):
                    if "libgfortran=" in d:
                        data["dependencies"][i] = "libgfortran"
                for i, p in enumerate(data["dependencies"][-1]["pip"]):
                    if "ruamel" in p:
                        data["dependencies"][-1]["pip"][i] = None
                v = [x for x in data["dependencies"][-1]["pip"] if x is not None]
                data["dependencies"][-1]["pip"] = v
            with open(yml_file, "w") as f:
                yaml.safe_dump(data, f)
        except:
            return

    def _bundle_environment_yml_has_ersilia(self, model_id: str) -> bool:
        search = "ersilia="
        dir = self._get_bundle_location(model_id)
        yml_file = os.path.join(dir, ENVIRONMENT_YML)
        if not os.path.exists(yml_file):
            return None
        with open(yml_file, "r") as f:
            data = yaml.safe_load(f)
            if not data["dependencies"]:
                return False
            for i, d in enumerate(data["dependencies"][:-1]):
                if search in d:
                    return True
            for i, p in enumerate(data["dependencies"][-1]["pip"]):
                if search in p:
                    return True
        return False

    def _bundle_dockerfile_has_ersilia(self, model_id: str) -> bool:
        dir = self._get_bundle_location(model_id)
        dockerfile = os.path.join(dir, DOCKERFILE)
        if not os.path.exists(dockerfile):
            return None
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "grep.txt")
        cmd = "grep -R 'ersilia' {0} > {1}".format(dockerfile, tmp_file)
        run_command(cmd)
        with open(tmp_file, "r") as f:
            grep = f.read()
        if grep:
            return True
        else:
            return False

    def _modify_bundle_dockerfile(self, model_id: str):
        if not self._bundle_uses_ersilia(model_id):
            return
        if self._bundle_environment_yml_has_ersilia(model_id):
            return
        if self._bundle_dockerfile_has_ersilia(model_id):
            return
        dockerfile = os.path.join(self._get_bundle_location(model_id), DOCKERFILE)
        text = ["", "# Install ersilia"]
        text += [
            "RUN pip install git+https://github.com/{0}/{1}.git".format(
                self.cfg.HUB.ORG, self.cfg.HUB.PACKAGE
            )
        ]  # TODO: add version with the @ character
        with open(dockerfile, "r") as f:
            lines = []
            for l in f:
                lines += [l.rstrip()]
        if lines[1][:10] == "MAINTAINER":
            lines = lines[:2] + text + lines[2:]
        else:
            lines = lines[:1] + text + lines[1:]
        with open(dockerfile, "w") as f:
            for l in lines:
                f.write(l + os.linesep)

    def _add_model_install_commands_to_requirements_txt(self, model_id: str):
        BundleRequirementsFile(model_id).add_model_install_commands()

    def _add_model_install_commands_to_environment_yml(self, model_id: str):
        BundleEnvironmentFile(model_id).add_model_install_commands()

    def _explicit_conda_python_path_in_run(self, model_id: str):
        dir = self._get_bundle_location(model_id)
        framework_dir = os.path.join(dir, model_id, "artifacts", "framework")
        if not os.path.exists(framework_dir):
            return
        run_files = []
        for l in os.listdir(framework_dir):
            if l.startswith("run") and l.endswith(".sh"):
                run_files += [l]
        if len(run_files) != 1:
            return
        run_file = os.path.join(framework_dir, run_files[0])
        self.logger.debug("Run file found in framework: {0}".format(run_file))
        with open(run_file, "r") as f:
            for l in f:
                if l.startswith("conda activate"):
                    self.logger.debug(
                        "A conda activate statement has been found. It is not advised to modify the conda path in this bash file."
                    )
                    return
        python_exec = SimpleConda().get_python_path_env(model_id)
        self.logger.debug("Python executable: {0}".format(python_exec))
        R = []
        with open(run_file, "r") as f:
            for r in f:
                if r.startswith("python "):
                    r = python_exec + r[6:]
                R += [r]
        with open(run_file, "w") as f:
            for r in R:
                f.write(r)

    def modify(self):
        """
        Modifies the model bundle to ensure compatibility and consistency.
        """
        # Add installs to requirements and environment
        self._add_model_install_commands_to_requirements_txt(self.model_id)
        self._add_model_install_commands_to_environment_yml(self.model_id)
        # Slightly modify bundle environment YAML file, if exists
        self._modify_bundle_environment_yml(self.model_id)
        # Slightly modify bundle Dockerfile, if necessary.
        self._modify_bundle_dockerfile(self.model_id)
        # Check if conda is really necessary, if not, use slim base docker image
        with open(
            os.path.join(self._model_path(self.model_id), PACKMODE_FILE), "r"
        ) as f:
            pack_mode = f.read()
        if pack_mode == "conda":
            needs_conda = True
            # Redirect python path from run command to conda path (if necessary)
            self._explicit_conda_python_path_in_run(self.model_id)
        else:
            needs_conda = BundleEnvironmentFile(self.model_id).needs_conda()
        dockerfile = BundleDockerfileFile(self.model_id)
        if needs_conda:
            self.logger.debug("Conda is needed")
            dockerfile.set_to_full()
        else:
            self.logger.debug("Conda is not needed")
            dockerfile.set_to_slim()
