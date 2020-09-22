"""Download utilities."""

import os
import zipfile
import requests
import tempfile
import base64
import shutil
from github import Github, GithubException


class StandardDownloader(object):

    def __init__(self):
        pass

    def _fetch(self, url, destination):
        pass

    def fetch_zip(self, url, destination):
        """A .zip file is assumed."""
        tmp_zip = tempfile.NamedTemporaryFile(dir=destination).name+".zip"
        self.download_file_from_google_drive(file_id, tmp_zip)
        with zipfile.ZipFile(tmp_zip, "r") as zip_ref:
            zip_ref.extractall(destination)
        os.remove(tmp_zip)

    def fetch_tarball(self, url, destination):
        """A tar.gz file is assumed"""

    def fetch(self, url, destination):
        """"""


class GoogleDriveDownloader(object):

    def __init__(self):
        pass

    @staticmethod
    def get_confirm_token(response):
        for key, value in response.cookies.items():
            if key.startswith('download_warning'):
                return value
        return None

    @staticmethod
    def save_response_content(response, destination):
        chunk_size = 32768
        with open(destination, "wb") as f:
            for chunk in response.iter_content(chunk_size):
                if chunk:
                    f.write(chunk)

    def download_file_from_google_drive(self, file_id, destination):
        url = "https://docs.google.com/uc?export=download"
        session = requests.Session()
        response = session.get(url, params={'id': file_id}, stream=True)
        token = self.get_confirm_token(response)
        if token:
            params = {'id': id, 'confirm': token}
            response = session.get(url, params=params, stream=True)
        self.save_response_content(response, destination)

    def fetch_zip(self, file_id, destination):
        """Download file from google docs. The file id is necessary. A .zip file is assumed."""
        tmp_zip = tempfile.NamedTemporaryFile(dir=destination).name+".zip"
        self.download_file_from_google_drive(file_id, tmp_zip)
        with zipfile.ZipFile(tmp_zip, "r") as zip_ref:
            zip_ref.extractall(destination)
        os.remove(tmp_zip)


class GitHubDownloader(object):

    def __init__(self, token=None):
        self.token = token
        if token is None:
            self.github = Github()
        else:
            self.github = Github(token)

    @staticmethod
    def get_sha_for_tag(repository, tag):
        """
        Returns a commit PyGithub object for the specified repository and tag.
        """
        branches = repository.get_branches()
        matched_branches = [match for match in branches if match.name == tag]
        if matched_branches:
            return matched_branches[0].commit.sha

        tags = repository.get_tags()
        matched_tags = [match for match in tags if match.name == tag]
        if not matched_tags:
            raise ValueError('No Tag or Branch exists with that name')
        return matched_tags[0].commit.sha

    def download_directory(self, repository, sha, server_path):
        """
        Download all contents at server_path with commit tag sha in
        the repository.
        """
        if os.path.exists(server_path):
            shutil.rmtree(server_path)

        os.makedirs(server_path)
        contents = repository.get_dir_contents(server_path, ref=sha)

        for content in contents:
            if content.type == 'dir':
                os.makedirs(content.path)
                self.download_directory(repository, sha, content.path)
            else:
                try:
                    path = content.path
                    file_content = repository.get_contents(path, ref=sha)
                    file_data = base64.b64decode(file_content.content)
                    file_out = open(content.path, "w+")
                    file_out.write(file_data)
                    file_out.close()
                except (GithubException, IOError) as exc:
                    print('Error processing %s: %s', content.path, exc)

    def fetch(self, folder, org, repo, tag, destination):
        organization = self.github.get_organization(org)
        repository = organization.get_repo(repo)
        sha = self.get_sha_for_tag(repository, tag)
        self.download_directory(repository, sha, folder)
