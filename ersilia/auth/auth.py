import os
from pathlib import Path

import yaml

try:
    from github import Github
    from github.GithubException import UnknownObjectException
except (ModuleNotFoundError, ImportError):
    Github = None
    UnknownObjectException = None


HOSTNAME = "github.com"
SECRET_REPO = "ersilia-os/ersilia-secrets"


class Auth:
    """
    This class handles authentication.
    """

    def __init__(self):
        self.hosts_yml = os.path.join(str(Path.home()), ".config", "gh", "hosts.yml")
        if os.path.exists(self.hosts_yml):
            with open(self.hosts_yml, "r") as f:
                self.hosts = yaml.safe_load(f)
        else:
            self.hosts = None
        self.hostname = HOSTNAME
        self.secret_repo = SECRET_REPO

    def login(self):
        """Login using GitHub"""
        pass
        #  TODO

    def logout(self):
        """Logout from GitHub"""
        pass
        # TODO

    def status(self):
        """See login status"""
        if self.hosts is None:
            return None
        else:
            return self.hosts[self.hostname]

    def user(self):
        """Get user"""
        if self.hosts is None:
            return None
        else:
            return self.hosts[self.hostname]["user"]

    def oauth_token(self):
        """Get OAuth Token"""
        if self.hosts is None:
            return None
        else:
            return self.hosts[self.hostname]["oauth_token"]

    def is_contributor(self):
        """Check if the GitHub user is a contributor of ersilia-os"""
        if Github is None:
            return False
        else:
            gh = Github(login_or_token=self.oauth_token())
            try:
                repo = gh.get_repo(self.secret_repo)
                if repo is None:
                    return False
                else:
                    return True
            except UnknownObjectException:
                return False
