import json
import logging
import os
import time
import sys

from github import Github, GithubException


class UpdateMetadata:
    """
    Class for reading the metadata file from a repo and updating it from the new model submission request
    """

    def __init__(self, issue_creator=None):
        self.log = self.logger()
        self.metadata_filename = "metadata.json"
        self.retries = 9
        self.retry_delay = 3
        self.token = os.environ.get("GITHUB_TOKEN")
        self.owner = os.environ.get("OWNER")
        self.repo = os.environ.get("REPO")
        self.branch = os.environ.get("BRANCH", "main")
        self.issue_creator = issue_creator
        self.metadata = None
        self.json_input = self.load_json_input()
        self.github = Github(self.token)

    def logger(self):
        """
        Create a logger for use in this class
        """
        self.log = logging.getLogger(__name__)
        self.log.setLevel(logging.INFO)
        formatter = logging.Formatter("%(levelname)s - %(message)s")
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        self.log.addHandler(handler)
        return self.log

    def load_json_input(self):
        """
        Helper function for loading the JSON input from the env vars
        This is the JSON data parsed from the new model submission request (GitHub issue)
        :return: dict
        """
        self.log.info(f"loading JSON input from env vars")

        # Load the JSON input from the env vars and convert it to a dict
        return json.loads(os.environ.get("JSON"))

    def read_metadata(self):
        """
        Read the metadata file from the repo
        :note: this is the metadata.json file in the repo that was just created by the new model submission request...
        ...and it should be empty (contain default values)
        """
        self.log.info(f"loading {self.metadata_filename} from {self.owner}/{self.repo}")

        # check to see if the repository has any contents
        for i in range(self.retries):
            try:
                # Get the repo and contents
                repo = self.github.get_repo(f"{self.owner}/{self.repo}")
                # fetch the repo contents
                contents = repo.get_contents(self.metadata_filename)
                break
            except GithubException as exception:
                # if the repo is empty, we will get a 404 error
                if exception.status != 404:
                    # if the error is not a 404, we will raise the exception
                    self.log.info(f"exception: {exception}")
                    raise exception

                # if we have reached the retry limit, we will raise the exception
                if i == self.retries:
                    self.log.info(f"retry limit reached ({self.retries}), exiting")
                    raise exception

                # if we have not reached the retry limit, we will wait and try again
                self.log.info(
                    f"waiting for {self.metadata_filename} to be created and trying again..."
                )
                time.sleep(self.retry_delay)

        # Load the metadata.json from the repo and convert it to a dict
        self.metadata = json.loads(contents.decoded_content.decode())

    def populate_metadata(self):
        """
        Populate the metadata file with the new model submission request
        We will mash the metadata file with the JSON input from the new model submission request
        """
        self.log.info(
            f"populating {self.metadata_filename} with new model submission request data"
        )

        # Match the metadata keys to the JSON input keys
        # We will only populate the metadata file with the JSON input if the metadata file is empty ("" or [])
        # This is gross, sorry
        if self.metadata["Identifier"] == "":
            self.metadata["Identifier"] = self.repo
        if self.metadata["Slug"] == "":
            self.metadata["Slug"] = self.json_input["slug"]
        if self.metadata["Title"] == "":
            self.metadata["Title"] = self.json_input["model_name"]
        if self.metadata["Description"] == "":
            self.metadata["Description"] = self.json_input["model_description"]
        if self.metadata["Publication"] == "":
            self.metadata["Publication"] = self.json_input["publication"]
        if self.metadata["Source Code"] == "":
            self.metadata["Source Code"] = self.json_input["source_code"]
        if self.metadata["License"] == "":
            self.metadata["License"] = self.json_input["license"]
        if self.metadata["Tag"] == []:
            # split the tags into a list andremove any whitespace
            tags = [tag.strip() for tag in self.json_input["tag"].split(",")]
            self.metadata["Tag"] = tags
        if self.metadata["Contributor"] == "":
            self.metadata["Contributor"] == self.issue_creator

    def write_metadata(self):
        """
        Write the metadata file to the repo
        In this function, we use the GitHub Python SDK to write the metadata file to the repo
        We will commit the file to the repo and push it to the main branch
        """
        self.log.info(
            f"committing {self.metadata_filename} to {self.owner}/{self.repo}"
        )

        # Get the repo and contents
        repo = self.github.get_repo(f"{self.owner}/{self.repo}")
        contents = repo.get_contents(self.metadata_filename)

        # Convert the metadata to a JSON formatted string
        metadata_string = json.dumps(self.metadata, indent=4).encode("utf-8")

        # Write the metadata to the repo and commit
        repo.update_file(
            contents.path,
            "Initialize metadata [skip ci]",
            metadata_string,
            contents.sha,
            branch=self.branch,
        )
        self.log.info(
            f"successfully wrote {self.metadata_filename} to {self.owner}/{self.repo}"
        )

    def run(self):
        """
        Run the script
        """
        self.read_metadata()
        self.populate_metadata()
        self.write_metadata()


if __name__ == "__main__":
    issue_creator = sys.argv[1]
    UpdateMetadata(issue_creator=issue_creator).run()
