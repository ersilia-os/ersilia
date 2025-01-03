import os
import shutil

from .. import ErsiliaBase
from ..default import GITHUB_ORG
from ..utils.dvc import DVCSetup
from ..utils.terminal import run_command
from . import EOS_TEMPLATE_REPOSITORY
from .rebase import TemplateRebaser


class ModelPublisher(ErsiliaBase):
    """
    Class for publishing models to GitHub.

    Parameters
    ----------
    model_id : str
        The ID of the model to be published.
    config_json : str
        Path to the configuration JSON file.
    credentials_json : str
        Path to the credentials JSON file.
    """

    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.model_id = model_id
        self.template_repo = EOS_TEMPLATE_REPOSITORY
        self.repo_remote_path = os.path.join(GITHUB_ORG, self.model_id)
        self.template_path = os.path.join(GITHUB_ORG, self.template_repo)
        self.cwd = os.path.abspath(os.getcwd())
        self.repo_path = os.path.join(self.cwd, self.model_id)
        self.message = "Commit from Ersilia"

    def rebase(self):
        """
        Rebase the model repository with the template repository.
        """
        rb = TemplateRebaser(
            model_id=self.model_id,
            template_repo=self.template_repo,
            config_json=self.config_json,
            credentials_json=self.credentials_json,
        )
        rb.rebase()
        shutil.copytree(rb.model_path, self.repo_path)
        self.logger.debug(
            "Rebased repository {0} with template {1}. Local path: {2}".format(
                self.model_id, self.template_repo, self.repo_path
            )
        )
        rb.clean()
        self.message = "Rebased from template {0}".format(self.template_repo)

    def create(self, public=True):
        """
        Create a new GitHub repository for the model.

        Parameters
        ----------
        public : bool, optional
            Whether the repository should be public or private. Default is True.
        """
        if public:
            flag = "--public"
        else:
            flag = "--private"
        cmd = "gh repo create {0} --confirm --template {1} {2}".format(
            self.repo_remote_path, self.template_path, flag
        )
        self.logger.debug(cmd)
        run_command(cmd)
        shutil.rmtree(self.repo_path)
        cmd = "gh repo clone {0}".format(self.repo_remote_path)
        run_command(cmd)
        self.logger.debug(
            "Created repository in GitHub. Local path: {0}".format(self.repo_path)
        )
        self.message = "Initial commit"

    def dvc(self):
        """
        Set up DVC (Data Version Control) for the model repository.
        """
        dvc = DVCSetup(local_repo_path=self.repo_path, model_id=self.model_id)
        dvc.gdrive_setup()
        dvc.set_dvc_gdrive()
        dvc.git_add_and_commit()

    def git_push(self, message=None):
        """
        Push changes to the GitHub repository.

        Parameters
        ----------
        message : str, optional
            The commit message. If not provided, a default message is used.
        """
        if not message:
            message = self.message
        os.chdir(self.repo_path)
        run_command("git add .")
        run_command("git commit -m '{0}'".format(message))
        run_command("git pull --rebase")
        run_command("git push")
        os.chdir(self.cwd)

    def push(self):
        """
        Set up DVC and push changes to the GitHub repository.
        """
        self.dvc()
        self.git_push()

    def test(self):
        """
        Test the publishing process.
        """
        pass

    def docker(self):
        """
        Handle Docker-related tasks.
        """
        pass
