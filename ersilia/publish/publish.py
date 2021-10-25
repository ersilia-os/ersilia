import os
import shutil

from ..utils.terminal import run_command

from .. import ErsiliaBase
from ..utils.dvc import DVCSetup
from . import EOS_TEMPLATE_REPOSITORY


class ModelPublisher(ErsiliaBase):
    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(
            self,
            config_json=config_json,
            credentials_json=credentials_json,
        )
        self.model_id = model_id
        self.repo_remote_path = os.path.join("ersilia-os", model_id)
        self.template_path = os.path.join("ersilia-os", EOS_TEMPLATE_REPOSITORY)
        self.cwd = os.path.abspath(os.getcwd())
        self.repo_path = os.path.join(self.cwd, self.model_id)

    def create(self, public=True):
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

    def dvc(self):
        dvc = DVCSetup(local_repo_path=self.repo_path, model_id=self.model_id)
        dvc.gdrive_setup()
        dvc.set_dvc_gdrive()
        dvc.git_add_and_commit()

    def git_push(self, message="Initial commit"):
        os.chdir(self.repo_path)
        run_command("git add .")
        run_command("git commit -m '{0}'".format(message))
        run_command("git pull --rebase")
        run_command("git push")
        os.chdir(self.cwd)

    def push(self):
        self.dvc()
        self.git_push()

    def test(self):
        pass

    def docker(self):
        pass
