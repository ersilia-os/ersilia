"""Download utilities"""

import os
import shutil
import subprocess
import tempfile
import uuid
import zipfile
from pathlib import Path

import requests

from .. import logger
from ..default import S3_BUCKET_URL
from ..utils.logging import make_temp_dir
from .terminal import run_command


class PseudoDownloader(object):
    """
    A class to simulate downloading by copying local directories.

    Parameters
    ----------
    overwrite : bool
        Whether to overwrite existing files.
    """

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
        """
        Copy entire directory.

        Parameters
        ----------
        src : str
            The source directory.
        dst : str
            The destination directory.
        """
        self._copy_local_directory(
            src, dst
        )  # TODO: Add smart functions to deal with zipped files etc.


class OsfDownloader(object):
    """
    A class to download files from the Open Science Framework (OSF).

    Parameters
    ----------
    overwrite : bool
        Whether to overwrite existing files.
    """

    def __init__(self, overwrite):
        self.overwrite = overwrite

    def fetch(self, project_id, filename, destination, tmp_folder):
        """
        Fetch a file from OSF.

        Parameters
        ----------
        project_id : str
            The OSF project ID.
        filename : str
            The name of the file to fetch.
        destination : str
            The destination directory.
        tmp_folder : str
            The temporary folder to use.
        """
        src = os.path.basename(filename)
        outfile = os.path.join(destination, src)
        if os.path.exists(outfile):
            if self.overwrite:
                os.remove(outfile)
            else:
                return
        cwd = os.getcwd()
        os.chdir(tmp_folder)
        run_command("osf -p %s fetch %s" % (project_id, filename))
        shutil.move(src, outfile)
        os.chdir(cwd)


class GoogleDriveDownloader(object):
    """
    A class to download files from Google Drive.
    """

    def __init__(self):
        pass

    @staticmethod
    def get_confirm_token(response):
        """
        Get the confirmation token from the response.

        Parameters
        ----------
        response : object
            The response object.

        Returns
        -------
        str
            The confirmation token, if available.
        """
        for key, value in response.cookies.items():
            if key.startswith("download_warning"):
                return value
        return None

    @staticmethod
    def save_response_content(response, destination):
        """
        Save the response content to a file.

        Parameters
        ----------
        response : object
            The response object.
        destination : str
            The destination file path.
        """
        chunk_size = 32768
        with open(destination, "wb") as f:
            for chunk in response.iter_content(chunk_size):
                if chunk:
                    f.write(chunk)

    def download_file_from_google_drive(self, file_id, destination):
        """
        Download a file from Google Drive.

        Parameters
        ----------
        file_id : str
            The Google Drive file ID.
        destination : str
            The destination path.
        """
        url = "https://docs.google.com/uc?export=download"
        session = requests.Session()
        response = session.get(url, params={"id": file_id}, stream=True)
        token = self.get_confirm_token(response)
        if token:
            params = {"id": id, "confirm": token}
            response = session.get(url, params=params, stream=True)
        self.save_response_content(response, destination)

    def fetch_zip(self, file_id, destination):
        """
        Download a ZIP file from Google Drive and extract it.

        Parameters
        ----------
        file_id : str
            The Google Drive file ID.
        destination : str
            The destination directory.
        """
        tmp_zip = tempfile.NamedTemporaryFile(dir=destination).name + ".zip"
        self.download_file_from_google_drive(file_id, tmp_zip)
        with zipfile.ZipFile(tmp_zip, "r") as zip_ref:
            zip_ref.extractall(destination)
        os.remove(tmp_zip)


class GitHubDownloader(object):
    """
    A class to download files and repositories from GitHub.

    Parameters
    ----------
    overwrite : bool
        Whether to overwrite existing files.
    token : str, optional
        The GitHub token for authentication. Default is None.
    """

    def __init__(self, overwrite, token=None):
        self.logger = logger

        self.overwrite = overwrite
        self.token = token

    @staticmethod
    def _repo_url(org, repo):
        return "https://github.com/{0}/{1}.git".format(org, repo)

    @staticmethod
    def _ungit(path):
        dotgit = os.path.join(path, ".git")
        if os.path.exists(dotgit):
            shutil.rmtree(dotgit)
        gitignore = os.path.join(path, ".gitignore")
        if os.path.exists(gitignore):
            os.remove(gitignore)
        gitattributes = os.path.join(path, ".gitattributes")
        if os.path.exists(gitattributes):
            os.remove(gitattributes)

    def _exists(self, destination):
        if os.path.exists(destination):
            return True
        else:
            return False

    def _clone_with_git(self, org, repo, destination):
        tmp_folder = os.path.abspath(make_temp_dir(prefix="ersilia-"))
        script = """
        #Disabling automatic LFS clone (s3 feature) by adding GIT_LFS_SKIP_SMUDGE=1
        cd {0}
        GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/{1}/{2}.git
        mv {2} {3}
        rm -r {0}
        """.format(tmp_folder, org, repo, destination)
        run_file = os.path.join(
            os.path.abspath(make_temp_dir(prefix="ersilia")), "run.sh"
        )
        with open(run_file, "w") as f:
            f.write(script)
        run_command("bash {0}".format(run_file))
        return self._exists(destination)

    def _list_lfs_files(self, destination):
        clean_lfs_files_list = []
        lfs_files = (
            subprocess.check_output(
                "cd " + destination + "; git lfs ls-files -l", shell=True
            )
            .decode("utf-8")
            .splitlines()
        )
        for lfs_file in lfs_files:
            # lfs_file takes a form of:
            #'7ca8461207e76ded1224f393e7bdb21973b5c6caece2a4e550e6396efe2cf9f7 - model/framework/RAscore/models/XGB_chembl_ecfp_counts/model.pkl'
            # where " - " separating SHA256 from a filename might also take a form of " * " if the file is already fetched from LFS.
            # In the next step we replace all the first occurences of " * " with " - " and do the split between the filename and sha256 of the file.
            file = lfs_file.replace(" * ", " - ", 1).split(" - ", 1)[1]
            sha256 = lfs_file.replace(" * ", " - ", 1).split(" - ", 1)[0]
            dic = {"file": file, "sha256": sha256}
            clean_lfs_files_list.append(dic)
        return clean_lfs_files_list

    def _download_s3_files(self, filename, repo, destination):
        file_url = f"{S3_BUCKET_URL}/{repo}/{filename}"
        local_filename = Path(destination) / filename

        try:
            response = requests.get(file_url, stream=True)
            response.raise_for_status()

            total_length = response.headers.get("content-length")
            total_length = int(total_length) if total_length else None

            local_filename.parent.mkdir(parents=True, exist_ok=True)

            with open(local_filename, "wb") as file:
                if total_length is None:
                    # No content-length header, write content directly
                    file.write(response.content)
                else:
                    self._download_large_file(response, file, total_length, filename)

            self.logger.info(
                f"✅ Successfully downloaded {filename} to {local_filename}"
            )
        except requests.RequestException as e:
            self.logger.error(
                f"❗ Could not download file {filename} from S3 bucket: {file_url}. "
                "Falling back to Git LFS if available."
            )
            self.logger.debug(f"Error details: {e}")
        except Exception as e:
            self.logger.error(f"❗ Unexpected error while downloading {filename}: {e}")

    def _download_large_file(self, response, file, total_length, filename):
        self.logger.info(f"Downloading large file {filename} from S3 bucket.")
        downloaded = 0

        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)
            downloaded += len(chunk)

    def _check_large_file_checksum(self, filename, destination):
        # This function takes filenames and checksums from lfs ls-files
        # and runs shasum -a 256 on actual files for comparison.
        expected_sha256 = filename["sha256"]
        actual_sha256 = (
            subprocess.check_output(
                "cd " + destination + "; shasum -a 256 " + filename["file"] + ";",
                shell=True,
            )
            .decode("utf-8")
            .split()[0]
        )

        if actual_sha256 == expected_sha256:
            return None
        else:
            self.logger.error(
                "❌ Checksum discrepancy in file {0}: expected {1}, actual {2}".format(
                    filename["file"], expected_sha256, actual_sha256
                )
            )
            return filename

    def _git_lfs(self, destination, filename):
        # This function is run for all the files that were
        # not downloaded from an S3 bucket or have unexpected sha256 value.
        self.logger.debug("⏳ Trying LFS clone for file {0}".format(filename))
        script = "cd {0}; git lfs pull --include {1}".format(destination, filename)
        tmp_folder = make_temp_dir(prefix="ersilia-")
        run_file = os.path.join(tmp_folder, "run_lfs.sh")
        with open(run_file, "w") as f:
            f.write(script)
        run_command("bash {0}".format(run_file))
        self.logger.success("✅")

    def _download_large_files(self, repo, destination):
        # This function downloads large files (as listed in .gitattributes).
        # T:
        # 1. gets the list of files to download from git lfs ls-files command.
        # 2. Tries to download listed files from S3 bucket.
        # 3. Checks sha256 for downloaded and existing files.
        # 4. If there's a discrepancy, tries to fetch only files in question
        # from git lfs.

        # Ad 1. List all LFS files.
        lfs_files_list = self._list_lfs_files(destination)
        files_with_incorrect_shasum = []
        # Download large files
        for filename in lfs_files_list:
            # Ad 2. Download file after file from lfs_files_list (list of dictionaries)
            self._download_s3_files(filename["file"], repo, destination)
            # Ad 3. Checks sha256 for downloaded and existing files.
            sha_check_result = self._check_large_file_checksum(filename, destination)
            # Ad 3. Add files with unexpected checksum to sha_check_result list.
            if sha_check_result is not None:
                files_with_incorrect_shasum.append(sha_check_result)
        # Ad 4. If there's a discrepancy, tries to fetch only conflicting files from LFS
        if files_with_incorrect_shasum:
            for files in files_with_incorrect_shasum:
                self._git_lfs(destination, files["file"])
        self.logger.info("🚀 Model starting...")

    def clone(self, org, repo, destination, ungit=False):
        """
        Clone a GitHub repository.

        Parameters
        ----------
        org : str
            The GitHub organization name.
        repo : str
            The GitHub repository name.
        destination : str
            The destination directory.
        ungit : bool, optional
            Whether to remove Git-related files. Default is False.
        """
        if os.path.exists(destination):
            if self.overwrite:
                shutil.rmtree(destination)
            else:
                return
        is_done = self._clone_with_git(org, repo, destination)
        if not is_done:
            raise Exception("Download from {0}/{1} did not work".format(org, repo))

        self._download_large_files(repo, destination)

        if ungit:
            self._ungit(destination)

    def _download_single_raw_by_branch(self, org, repo, branch, repo_path, destination):
        url = "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
            org, repo, branch, repo_path
        )
        response = requests.get(url)
        if response.status_code != 200:
            return False
        with open(destination, "wb") as f:
            f.write(response.content)
        return True

    def _download_single_raw(self, org, repo, repo_path, destination):
        is_done = self._download_single_raw_by_branch(
            org, repo, "master", repo_path, destination
        )
        if is_done:
            return True
        else:
            is_done = self._download_single_raw_by_branch(
                org, repo, "main", repo_path, destination
            )
        return is_done

    def download_single(self, org, repo, repo_path, destination):
        """
        Download a single file from a GitHub repository.

        Parameters
        ----------
        org : str
            The GitHub organization name.
        repo : str
            The GitHub repository name.
        repo_path : str
            The path to the file in the repository.
        destination : str
            The destination path.

        Returns
        -------
        bool
            True if the file was downloaded successfully, False otherwise.
        """
        if os.path.exists(destination):
            if self.overwrite:
                if os.path.isfile(destination):
                    os.remove(destination)
                if os.path.isdir(destination):
                    shutil.rmtree(destination)
            else:
                return
        is_done = self._download_single_raw(org, repo, repo_path, destination)
        if is_done:
            return True
        else:
            tmpdir = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
            self.clone(org, repo, tmpdir)
            source = os.path.join(tmpdir, repo_path)
            if os.path.exists(source):
                if os.path.isfile(source):
                    shutil.copyfile(source, destination)
                if os.path.isdir(source):
                    shutil.copytree(source, destination)
            shutil.rmtree(tmpdir)
        return self._exists(destination)


class S3Downloader(object):
    """
    A class to download files from an S3 bucket.
    """

    def __init__(self):
        pass

    def download_from_s3(self, bucket_url, file_name, destination):
        """
        Download a file from an S3 bucket.

        Parameters
        ----------
        bucket_url : str
            The URL of the S3 bucket.
        file_name : str
            The name of the file to download.
        destination : str
            The destination path.
        """
        s3_url = bucket_url + "/" + file_name
        response = requests.get(s3_url)
        response.raise_for_status()
        with open(destination, "wb") as f:
            f.write(response.content)
