from .. import ErsiliaBase
import requests
import subprocess
import time
import json
from urllib.request import urlopen
from ..hub.content.card import RepoMetadataFile


class ModelInspector(ErsiliaBase):
    def __init__(self, model, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model = model

    def check_repo_exists(self, flag):
        """
        Verify that the repository exists at a given link.

        Args:
            flag (int): Flag indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the flag, returns boolean or a message indicating the check status.
        """
        url = f"https://github.com/ersilia-os/{self.model}"
        response = requests.head(url)
        if flag == 0:
            return response.status_code == 200
        elif response.status_code == 200:
            return "Check passed."
        return f"Connection invalid, no github repository found at https://github.com/ersilia-os/{self.model}. Please check that this repository exists in the ersilia-os database."

    def check_complete_metadata(self, flag):
        """
        Search for specific keys in metadata JSON file.

        Args:
            flag (int): Flag indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the flag, returns boolean or a message indicating the check status.
        """
        url = f"https://raw.githubusercontent.com/ersilia-os/{self.model}/main/metadata.json"  # Get raw file from GitHub
        if requests.head(url).status_code != 200:  # Make sure repo exists
            if flag == 0:
                return False
            return f"Metadata file could not be loated for model {self.model}. Please check that the link https://raw.githubusercontent.com/ersilia-os/{self.model}/main/metadata.json is valid."

        response = requests.get(url)
        file = response.json()  # Save as json object

        check_passed = True
        details = ""

        test = RepoMetadataFile(self.model)

        try:
            RepoMetadataFile.read_information(test)

        except Exception as e:
            details = (
                details + f"Error encountered when parsing through metadata file: {e} "
            )
            check_passed = False

        if file is not None:
            try:
                if (
                    file["Publication"]
                    and file["Source Code"]
                    and file["S3"]
                    and file["DockerHub"]
                ):  # Parse through json object and ensure

                    pub_url_works = (
                        requests.head(file["Publication"]).status_code != 404
                    )
                    if not pub_url_works:
                        details = (
                            details
                            + f"The URL ({file['Publication']}) listed in the publication field of the metadata was not found (Error code {requests.head(file['Publication']).status_code}). Please verify that this link is accurate. "
                        )
                        check_passed = False

                    source_url_works = (
                        requests.head(file["Source Code"]).status_code != 404
                    )
                    if not source_url_works:
                        details = (
                            details
                            + f"The URL ({file['Source Code']}) listed in the source code field of the metadata was not able to be accessed (Error code {requests.head(file['Source Code']).status_code}). Please verify that this link is accurate. "
                        )
                        check_passed = False

                    s3_url_works = requests.head(file["S3"]).status_code != 404
                    if not s3_url_works:
                        details = (
                            details
                            + f"The URL ({file['S3']}) listed in the S3 field of the metadata was not able to be accessed (Error code {requests.head(file['S3']).status_code}). Please verify that this link is accurate. "
                        )
                        check_passed = False

                    docker_url_works = (
                        requests.head(file["DockerHub"]).status_code != 404
                    )
                    if not docker_url_works:
                        details = (
                            details
                            + f"The URL ({file['DockerHub']}) listed in the DockerHub field of the metadata was not able to be accessed (Error code {requests.head(file['DockerHub']).status_code}). Please verify that this link is accurate. "
                        )
                        check_passed = False

                else:
                    details = (
                        details
                        + "A field in the metadata is empty. Please check that there are values listed for Publication, Source Code, S3, and DockerHub listed in the metadata. "
                    )
                    check_passed = False

            except KeyError:  # If a given key not present in json file return false
                check_passed = False
                details = (
                    details
                    + "Not all required fields were found in the metadata. Please check that Publication, Source Code, S3, and DockerHub are all present in the metadata. "
                )

            except requests.exceptions.ConnectionError:
                check_passed = False
                details = (
                    details
                    + "Connection failed when trying to access a URL listed in the metadata. Please check that URLs are all accurate and accessible. "
                )

        if flag == 0:
            return check_passed
        if details:
            return details
        return "Check passed."

    def check_complete_folder_structure(self, flag):
        """
        Validate folder structure of the repository.

        Args:
            flag (int): Flag indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the flag, returns boolean or a message indicating the check status.
        """
        check_passed = True
        details = ""
        url = f"https://github.com/ersilia-os/{self.model}"
        if requests.head(url).status_code != 200:  # Make sure repo exists
            if flag == 0:
                return False
            return f"Repository could not be loated for model {self.model}. Please check that the link https://github.com/ersilia-os/{self.model} is valid."

        folders = [
            ".github/workflows",
            "model",
            "src",
            "model/checkpoints",
            "model/framework",
        ]
        for name in folders:
            response = requests.get(
                url + "/tree/main/" + name
            )  # Check if the folders are present in a given repository
            if response.status_code != 200:
                details = (
                    details
                    + f"No {name} folder could be found. please check that the link {url}/blob/main/{name} is valid. "
                )
                check_passed = False  # If the folder URL is not valid return false

        files = ["LICENSE", "Dockerfile"]
        for name in files:
            response = requests.get(
                url + "/blob/main/" + name
            )  # Check if the files are present in a given repository
            if response.status_code != 200:
                details = (
                    details
                    + f"No {name} file could be found. please check that the link {url}/blob/main/{name} is valid. "
                )
                check_passed = False  # If the folder URL is not valid return false

        if flag == 0:
            return check_passed
        if check_passed:
            return "Check passed."
        return details

    def check_dependencies_are_valid(self, flag):
        """
        Check dependencies specified in the Dockerfile.

        Args:
            flag (int): Flag indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the flag, returns boolean or a message indicating the check status.
        """
        check_passed = True
        details = ""
        url = f"https://raw.githubusercontent.com/ersilia-os/{self.model}/main/Dockerfile"  # Get raw file from GitHub
        if requests.head(url).status_code != 200:  # Make sure repo exists
            if flag == 0:
                return False
            return f"Dockerfile could not be loated for model {self.model}. Please check that the link https://raw.githubusercontent.com/ersilia-os/{self.model}/main/Dockerfile is valid."

        response = requests.get(url)
        file = response.text
        lines = file.split("\n")
        lines = [s for s in lines if s]

        for line in lines:
            if line.startswith("RUN pip install"):
                info = line.split("==")
                install = line.split("RUN ")
                result = subprocess.run(
                    install,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                if result.returncode != 0 and flag == 0:
                    return False
                if result.returncode != 0 and flag == 1:
                    return f"Error running {install}, {result.stderr}"

                if len(info) < 2:
                    details = (
                        details
                        + f"No version specification found for the line {info[0]} in the Docker file. "
                    )
                    check_passed = False
                else:
                    specification = info[1]
                    if specification.strip() == "":
                        details = (
                            details
                            + f"No version specification found for the line {info[0]} in the Docker file. "
                        )
                        check_passed = False

        if "WORKDIR /repo" not in lines[len(lines) - 2]:
            details = (
                details + f"Dockerfile is missing 'WORKDIR /repo' in the right place. "
            )
            check_passed = False

        if (
            "COPY . /repo" not in lines[len(lines) - 1]
            and "COPY ./repo" not in lines[len(lines) - 1]
        ):
            details = (
                details
                + f"Dockerfile is missing 'COPY . /repo' in the right place or has incorrect syntax. "
            )
            check_passed = False

        if flag == 0:
            return check_passed
        if check_passed:
            return "Check passed."
        return details

    def check_comptuational_performance(self, verbose):
        """
        Measure computational performance by serving the model and running predictions.

        Args:
            verbose (int): Verbose indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the verbose, returns boolean or a message indicating the check status.
        """
        details = ""

        for n in (1, 10, 100):

            cmd = (
                f"ersilia serve {self.model} && "
                f"ersilia example -f my_input.csv -n {n} && "
                "ersilia run -i my_input.csv && "
                "ersilia close"
            )

            startTime = time.time()

            process = subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if process.returncode != 0:
                if verbose == 0:
                    return False
                return f"Error serving model: {process.stdout}, {process.stderr}"

            endTime = time.time()

            executionTime = endTime - startTime
            details += f"Execution time ({n} Prediction(s)): {executionTime} seconds. "

        if verbose == 1:
            return details
        return True

    def check_no_extra_files(self, flag):
        """
        Ensure that there are no excess files in the root directory of the repository.

        Args:
            flag (int): Flag indicating whether to return boolean or detailed message.

        Returns:
            bool or str: Depending on the flag, returns boolean or a message indicating the check status.
        """
        check_passed = True
        details = ""
        url = f"https://api.github.com/repos/ersilia-os/{self.model}/contents"
        if requests.head(url).status_code != 200:  # Make sure repo exists
            if flag == 0:
                return False
            return f"Repository could not be loated for model {self.model}. Please check that the link https://api.github.com/repos/ersilia-os/{self.model}/contents is valid."
        headers = {"Accept": "application/vnd.github.v3+json"}
        response = requests.get(url)

        folders = [
            ".github",
            "model",
            "src",
            ".gitignore",
            "Dockerfile",
            "LICENSE",
            "README.md",
            "metadata.json",
            "pack.py",
            ".gitattributes",
        ]
        for item in response.json():
            name = item["name"]
            if name not in folders:
                details = (
                    details
                    + f"Unexpected folder/file {name} found in root directory. Please check that {name} is a valid folder/file. "
                )
                check_passed = False  # If the folder URL is not valid return false

        if flag == 0:
            return check_passed
        if check_passed:
            return "Check passed."
        return details
