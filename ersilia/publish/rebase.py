import os
import shutil

from .. import ErsiliaBase
from ..default import GITHUB_ORG
from ..utils.terminal import run_command


class _FileFolderRebaser(object):
    def __init__(self, model_path, template_path):
        self.model_path = model_path
        self.template_path = template_path

    def do_file(self, name, overwrite):
        src = os.path.join(self.template_path, name)
        trg = os.path.join(self.model_path, name)
        if not os.path.exists(src):
            return
        if overwrite:
            shutil.copyfile(src, trg)
        else:
            if os.path.exists(trg):
                return
            else:
                shutil.copyfile(src, trg)

    def do_folder(self, name, overwrite):
        src = os.path.join(self.template_path, name)
        trg = os.path.join(self.model_path, name)
        if not os.path.exists(src):
            return
        if overwrite:
            shutil.copytree(src, trg)
        else:
            if os.path.exists(trg):
                return
            else:
                shutil.copytree(src, trg)


class TemplateRebaser(ErsiliaBase):
    """
    Class for rebasing model repositories with a template repository.

    Parameters
    ----------
    model_id : str
        The ID of the model to be rebased.
    template_repo : str, optional
        The name of the template repository. Default is 'eos-template'.
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(
        self,
        model_id: str,
        template_repo="eos-template",
        config_json=None,
        credentials_json=None,
    ):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
        self.template_repo = template_repo
        self.root = os.path.abspath(self._tmp_dir)
        self.cwd = os.getcwd()
        self.model_path = os.path.join(self.root, self.model_id)
        self.template_path = os.path.join(self.root, self.template_repo)
        self.clean()
        self.file_folder_rebaser = _FileFolderRebaser(
            self.model_path, self.template_path
        )

    def clone_template(self):
        """
        Clone the template repository.
        """
        os.chdir(self.root)
        run_command("gh repo clone {0}/{1}".format(GITHUB_ORG, self.template_repo))
        os.chdir(self.cwd)

    def clone_current_model(self):
        """
        Clone the current model repository.
        """
        os.chdir(self.root)
        run_command("gh repo clone {0}/{1}".format(GITHUB_ORG, self.model_id))
        os.chdir(self.cwd)

    def dvc_part(self):
        """
        Set up DVC (Data Version Control) for the model repository.
        """
        self.file_folder_rebaser.do_file("data.h5", overwrite=False)
        self.file_folder_rebaser.do_file(".dvcignore", overwrite=False)
        self.file_folder_rebaser.do_folder(".dvc", overwrite=False)

    def clean(self):
        """
        Clean up temporary directories.
        """
        if os.path.exists(self.template_path):
            self.logger.debug("Cleaning {0}".format(self.template_path))
            run_command("rm -rf {0}".format(self.template_path))
        if os.path.exists(self.model_path):
            self.logger.debug("Cleaning {0}".format(self.model_path))
            run_command("rm -rf {0}".format(self.model_path))

    # TODO: Add other rebasing options
    def rebase(self):
        """
        Rebase the model repository with the template repository.
        """
        self.clone_template()
        self.clone_current_model()
        self.dvc_part()
