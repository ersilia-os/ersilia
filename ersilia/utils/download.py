"""Download utilities"""

import os
import zipfile
import requests
import shutil
import tempfile
import uuid
import requests
import sys
import subprocess
from click import echo
from .terminal import run_command
from .. import logger

from ..default import S3_BUCKET_URL, S3_BUCKET_URL_ZIP


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
        self._copy_local_directory(
            src, dst
        )  # TODO: Add smart functions to deal with zipped files etc.


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
        run_command("osf -p %s fetch %s" % (project_id, filename))
        shutil.move(src, outfile)
        os.chdir(cwd)


class GoogleDriveDownloader(object):
    def __init__(self):
        pass

    @staticmethod
    def get_confirm_token(response):
        for key, value in response.cookies.items():
            if key.startswith("download_warning"):
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
        response = session.get(url, params={"id": file_id}, stream=True)
        token = self.get_confirm_token(response)
        if token:
            params = {"id": id, "confirm": token}
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
        tmp_folder = os.path.abspath(tempfile.mkdtemp(prefix="ersilia-"))
        script = """
        #Disabling automatic LFS clone (s3 feature) by adding GIT_LFS_SKIP_SMUDGE=1
        cd {0}
        GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/{1}/{2}.git
        mv {2} {3}
        rm {0}
        """.format(
            tmp_folder, org, repo, destination
        )
        run_file = os.path.join(
            os.path.abspath(tempfile.mkdtemp(prefix="ersilia")), "run.sh"
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
        # This function takes S3 filename as input and tries to download it
        # from a location given in S3_BUCKET_URL at default.py

        file_url = S3_BUCKET_URL + "/" + repo + "/" + filename
        local_filename = destination + "/" + filename
        try:
            with requests.get(file_url, stream=True) as r:
                r.raise_for_status()
                dl = 0
                total_length = int(r.headers.get("content-length"))
                if total_length is None:  # no content length header
                    f.write(r.content)
                else:
                    with open(local_filename, "wb") as f:
                        echo(
                            "Downloading large file {} from S3 bucket.".format(filename)
                        )
                        for chunk in r.iter_content(chunk_size=8192):
                            dl += len(chunk)
                            f.write(chunk)
                            done = int(50 * dl / total_length)
                            sys.stdout.write(
                                "\r[%s%s]" % ("=" * done, " " * (50 - done))
                            )
                            sys.stdout.flush()
                    echo("✅\n")
        except:
            self.logger.error(
                "❗Could not download file {} from S3 bucket.\n We will try Git LFS.".format(
                    file_url
                )
            )

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
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
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
    def __init__(self):
        pass

    def download_from_s3(self, bucket_url, file_name, destination):
        s3_url = bucket_url + "/" + file_name
        response = requests.get(s3_url)
        response.raise_for_status()
        with open(destination, "wb") as f:
            f.write(response.content)
