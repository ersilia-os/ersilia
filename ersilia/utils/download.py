"""Download utilities"""

import os
import zipfile
import requests
import shutil
import tempfile
import uuid
from .terminal import run_command


class PseudoDownloader(object):

    def __init__(self, overwrite):
        self.overwrite = overwrite

    def _fetch(self, url, destination):
        pass

    def _copy_local_directory(self, src, dst):
        if os.path.exists(dst):
            if self.overwrite:
                shutil.rmtree(dst)
            else:
                return
        shutil.copytree(src, dst)

    def fetch(self, src, dst):
        """Copy entire directory"""
        self._copy_local_directory(src, dst)  # TODO: Add smart functions to deal with zipped files etc.


class OsfDownloader(object):

    def __init__(self, overwrite):
        self.overwrite = overwrite

    def fetch(self, project_id, filename, destination, tmp_folder):
        src = os.path.basename(filename)
        outfile = os.path.join(destination, src)
        if os.path.exists(outfile):
            if self.overwrite:
                os.remove(outfile)
            else:
                return
        cwd = os.getcwd()
        os.chdir(tmp_folder)
        run_command("osf -p %s fetch %s" % (project_id, filename), quiet=True)
        shutil.move(src, outfile)
        os.chdir(cwd)


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
        tmp_zip = tempfile.NamedTemporaryFile(dir=destination).name + ".zip"
        self.download_file_from_google_drive(file_id, tmp_zip)
        with zipfile.ZipFile(tmp_zip, "r") as zip_ref:
            zip_ref.extractall(destination)
        os.remove(tmp_zip)


class GitHubDownloader(object):

    def __init__(self, overwrite, token=None):
        self.overwrite = overwrite
        self.token = token

    @staticmethod
    def _repo_url(org, repo):
        return "https://github.com/{0}/{1}.git".format(org, repo)

    def clone(self, org, repo, destination):
        if os.path.exists(destination):
            if self.overwrite:
                shutil.rmtree(destination)
            else:
                return
        if shutil.which("gh") is not None:
            cmd = "gh repo clone {0}/{1} {2}".format(org, repo, destination)
            run_command(cmd, quiet=True)
        else:
            import pygit2
            pygit2.clone_repository(url=self._repo_url(org, repo), path=destination)

    def download_single(self, org, repo, repo_path, destination):
        if os.path.exists(destination):
            if self.overwrite:
                if os.path.isfile(destination):
                    os.remove(destination)
                if os.path.isdir(destination):
                    shutil.rmtree(destination)
            else:
                return
        tmpdir = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        self.clone(org, repo, tmpdir)
        source = os.path.join(tmpdir, repo_path)
        if os.path.exists(source):
            if os.path.isfile(source):
                shutil.copyfile(source, destination)
            if os.path.isdir(source):
                shutil.copytree(source, destination)
        shutil.rmtree(tmpdir)
