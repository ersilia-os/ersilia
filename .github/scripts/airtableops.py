import os

import pyairtable
import requests
import yaml

from ersilia.hub.content.card import BaseInformation, RepoMetadataFile
from ersilia.utils.logging import make_temp_dir
from ersilia.utils.terminal import run_command

GITHUB_ORG = "ersilia-os"
AIRTABLE_MODEL_HUB_BASE_ID = "appgxpCzCDNyGjWc8"
AIRTABLE_MODEL_HUB_TABLE_NAME = "Models"
AIRTABLE_MAX_ROWS = 100000
AIRTABLE_PAGE_SIZE = 100
ROK_URL = "https://ersilia-model-hub.s3.eu-central-1.amazonaws.com/read_only_keys.json"


class AirtableInterface:
    def __init__(self, mode="ro", api_key=None):
        self.base_id = AIRTABLE_MODEL_HUB_BASE_ID
        self.table_name = AIRTABLE_MODEL_HUB_TABLE_NAME
        self.max_rows = AIRTABLE_MAX_ROWS
        self.page_size = AIRTABLE_PAGE_SIZE
        if mode == "ro":
            self.table = self._get_ro_table()
        elif mode == "rw":
            if api_key is None:
                raise ValueError("API key must be provided for read-write mode")
            self.table = self._get_rw_table(api_key)
        else:
            raise ValueError("Mode must be either 'ro' or 'rw'")

    def _create_table(self, api_key):
        return pyairtable.Table(api_key, self.base_id, self.table_name)

    @staticmethod
    def _get_ro_airtable_api_key():
        r = requests.get(ROK_URL)
        data = r.json()
        return data["AIRTABLE_READONLY_API_KEY"]

    def _get_rw_table(self, api_key):
        return self._create_table(api_key=api_key)

    def _get_ro_table(self):
        ro_api_key = self._get_ro_airtable_api_key()
        return self._create_table(api_key=ro_api_key)

    def items(self):
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                yield record

    def items_all(self):
        records = self.table.all()
        for record in records:
            yield record


class AirtableMetadata(AirtableInterface):
    def __init__(self, model_id, api_key=None, mode="ro"):
        super().__init__(mode=mode, api_key=api_key)
        self.model_id = model_id
        self._empty_row_message = "The AirTable field Identifier was not found! Please check that there are not empty rows."

    def _find_record(self):
        data = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                try:
                    if self.model_id == record["fields"]["Identifier"]:
                        data = record["fields"]
                except:
                    print(self._empty_row_message)

        return data

    def _find_airtable_record_id(self):
        rec_id = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                try:
                    if self.model_id == record["fields"]["Identifier"]:
                        rec_id = record["id"]
                except:
                    print(self._empty_row_message)
        return rec_id

    def read_information(self):
        data = self._find_record()
        bi = BaseInformation()
        bi.from_dict(data)
        return bi

    def write_information(self, data: BaseInformation):
        d = data.as_dict()
        d["GitHub"] = data.github
        rec_id = self._find_airtable_record_id()
        if rec_id is None:
            print(
                "Model {0} does not exist in AirTable. Creating new record".format(
                    self.model_id
                )
            )
            self.table.create(d)
        else:
            print("Model {0} exists in AirTable. Updating record".format(self.model_id))
            self.table.update(record_id=rec_id, fields=d)


class ReadmeMetadata:
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id

    def read_information(self):
        print("Cannot read directly from README file. Using AirTable instead")
        am = AirtableMetadata(model_id=self.model_id)
        bi = am.read_information()
        print(bi.as_dict())
        return bi

    def write_information(self, data: BaseInformation, readme_path=None):
        d = data.as_dict()
        text = "# {0}\n\n".format(d["Title"])
        text += "{0}\n\n".format(d["Description"].rstrip("\n"))
        text += "## Identifiers\n\n"
        text += "* EOS model ID: `{0}`\n".format(d["Identifier"])
        text += "* Slug: `{0}`\n\n".format(d["Slug"])
        text += "## Characteristics\n\n"
        text += "* Input: `{0}`\n".format(", ".join(d["Input"]))
        text += "* Input Shape: `{0}`\n".format(d["Input Shape"])
        text += "* Task: `{0}`\n".format(", ".join(d["Task"]))
        text += "* Output: `{0}`\n".format(", ".join(d["Output"]))
        text += "* Output Type: `{0}`\n".format(", ".join(d["Output Type"]))
        text += "* Output Shape: `{0}`\n".format(d["Output Shape"])
        text += "* Interpretation: {0}\n\n".format(d["Interpretation"])
        text += "## Baseline Performance\n\n"
        text += "* Computational Performance For One Input: `{0}`\n".format(d["Computational Performance 1"])
        text += "* Computational Performance For Ten Input: `{0}`\n".format(d["Computational Performance 10"])
        text += "* Computational Performance For Hundred Input: `{0}`\n".format(d["Computational Performance 100"])
        text += "## References\n\n"
        text += "* [Publication]({0})\n".format(d["Publication"])
        text += "* [Source Code]({0})\n".format(d["Source Code"])
        text += "* Ersilia contributor: [{0}](https://github.com/{0})\n\n".format(
            d["Contributor"]
        )
        text += "## Ersilia model URLs\n"
        text += "* [GitHub]({0})\n".format(data.github)
        if "S3" in d:
            text += "* [AWS S3]({0})\n".format(d["S3"])
        if "DockerHub" in d:
            text += "* [DockerHub]({0}) ({1})\n".format(
                d["DockerHub"], ", ".join(d["Docker Architecture"])
            )
        text += "\n"
        text += "## Citation\n\n"
        text += "If you use this model, please cite the [original authors]({0}) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).\n\n".format(
            d["Publication"]
        )
        text += "## License\n\n"
        text += "This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a {0} license.\n\n".format(
            d["License"]
        )
        text += "Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.\n\n"
        text += "## About Us\n\n"
        text += "The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.\n\n"
        text += "[Help us](https://www.ersilia.io/donate) achieve our mission!"
        if readme_path is None:
            return text
        else:
            with open(readme_path, "w") as f:
                f.write(text)


class FileUpdater:
    def __init__(self, model_id=None, repo_path=None, commit=True):
        self.model_id = model_id
        if repo_path is not None:
            self.repo_path = os.path.abspath(repo_path)
        else:
            self.repo_path = None
        self.commit = commit
        self.tmp_folder = make_temp_dir(prefix="ersilia-os")
        self.cwd = os.getcwd()

    def _git_clone(self):
        run_command(
            "cd {0}; GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/{1}/{2}; cd {3}".format(
                self.tmp_folder, GITHUB_ORG, self.model_id, self.cwd
            )
        )

    def _git_push(self):
        if self.repo_path is None:
            run_command(
                'cd {0}/{1}; GIT_LFS_SKIP_SMUDGE=1 git add .; git commit -m "Updating README file from AirTable [skip ci]"; git push; cd {2}'.format(
                    self.tmp_folder, self.model_id, self.cwd
                )
            )
        else:
            run_command(
                'cd {0}; GIT_LFS_SKIP_SMUDGE=1 git add .; git commit -m "Updating README file from AirTable [skip ci]"; git push; cd {1}'.format(
                    self.repo_path, self.cwd
                )
            )


class ReadmeUpdater(FileUpdater):
    def __init__(self, model_id, repo_path=None):
        super().__init__(model_id=model_id, repo_path=repo_path)

    def update_remote(self):
        self._git_clone()
        rm = ReadmeMetadata(model_id=self.model_id)
        bi = rm.read_information()
        tmp_file = os.path.join(self.tmp_folder, self.model_id, "README.md")
        rm.write_information(data=bi, readme_path=tmp_file)
        if self.commit:
            self._git_push()

    def update_local(self):
        rm = ReadmeMetadata(model_id=self.model_id)
        bi = rm.read_information()
        readme_file = os.path.join(self.repo_path, "README.md")
        rm.write_information(data=bi, readme_path=readme_file)
        if self.commit:
            self._git_push()

    def update(self):
        if self.repo_path is None:
            self.update_remote()
        else:
            self.update_local()


class MetadataUpdater(FileUpdater):
    def __init__(self, model_id, repo_path=None):
        super().__init__(model_id=model_id, repo_path=repo_path)

    def read_information(self):
        am = AirtableMetadata(model_id=self.model_id)
        bi = am.read_information()
        print(bi.as_dict())                                                                     
        return bi

    def update_remote(self):
        self._git_clone()
        am = AirtableMetadata(model_id=self.model_id)
        bi = am.read_information()
        rm = RepoMetadataFile(model_id=self.model_id, config_json=None)
        rm.get_json_or_yaml_file()
        if rm.is_json:
            tmp_file = os.path.join(self.tmp_folder, self.model_id, "metadata.json")
        else:
            tmp_file = os.path.join(self.tmp_folder, self.model_id, "metadata.yml")
        rm.write_information(data=bi, json_or_yaml_path=tmp_file)
        if self.commit:
            self._git_push()

    def update_local(self):
        am = AirtableMetadata(model_id=self.model_id)
        bi = am.read_information()
        rm = RepoMetadataFile(model_id=self.model_id, config_json=None)
        rm.get_json_or_yaml_file()
        if rm.is_json:
            tmp_file = os.path.join(self.tmp_folder, self.model_id, "metadata.json")
        else:
            tmp_file = os.path.join(self.tmp_folder, self.model_id, "metadata.yml")
        rm.write_information(data=bi, json_or_yaml_path=tmp_file)
        if self.commit:
            self._git_push()

    def update(self):
        if self.repo_path is None:
            self.update_remote()
        else:
            self.update_local()


def insert_metadata_to_airtable(model, contributor, api_key):
    # Works with airtable-insert option
    METADATA_YAML_FILE = "metadata.yml"
    ORG = "ersilia-os"
    BRANCH = "main"
    yaml_path = (
        f"https://raw.githubusercontent.com/{ORG}/{model}/{BRANCH}/{METADATA_YAML_FILE}"
    )
    github = f"https://github.com/ersilia-os/{model}"
    ai = AirtableInterface(mode="rw", api_key=api_key)

    r = requests.get(yaml_path)
    if r.status_code == 200:
        text = r.content
        data = yaml.safe_load(text)

    airtable_data = {}
    airtable_data["Identifier"] = model
    airtable_data["Slug"] = data["Slug"]
    airtable_data["Title"] = data["Title"]
    airtable_data["GitHub"] = github
    airtable_data["Contributor"] = contributor
    airtable_data["Status"] = "In progress"
    if data["Publication"] != "":
        airtable_data["Publication"] = data["Publication"]
    if data["Source Code"] != "":
        airtable_data["Source Code"] = data["Source Code"]
    ai.table.create(airtable_data)

def update_metadata_to_airtable(user, repo, branch, api_key):
    # Works with airtable-update option
    rm = RepoMetadataFile(model_id=repo, config_json=None)
    data = rm.read_information(org=user, branch=branch)
    am = AirtableMetadata(model_id=repo, api_key=api_key, mode="rw")
    am.write_information(data)

def update_readme_from_airtable(repo, path):
    # Works with readme-update option
    rm = ReadmeUpdater(model_id=repo, repo_path=path, commit=False)
    rm.update()

def update_metadata_from_airtable(repo, path):
    # Works with metadata-update option
    md = MetadataUpdater(model_id=repo, repo_path=path)
    md.update()

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()

    subparsers = parser.add_subparsers(dest="command")

    # Main commands
    airtable_insert = subparsers.add_parser(
        "airtable-insert", help="Insert metadata to AirTable"
    )
    airtable_update = subparsers.add_parser(
        "airtable-update", help="Update metadata to AirTable"
    )
    readme_update = subparsers.add_parser(
        "readme-update", help="Update README from AirTable"
    )
    metadata_update = subparsers.add_parser(
        "metadata-update", help="Update metadata from AirTable"
    )

    # Options for airtable-insert
    airtable_insert.add_argument("--model", type=str, required=True)
    airtable_insert.add_argument("--contributor", type=str, required=True)
    airtable_insert.add_argument("--api-key", type=str, required=True)

    # Options for airtable-update
    airtable_update.add_argument("--user", type=str, required=True)
    airtable_update.add_argument("--repo", type=str, required=True)
    airtable_update.add_argument("--branch", type=str, required=True)
    airtable_update.add_argument("--api-key", type=str, required=True)

    # Options for readme-update
    readme_update.add_argument("--repo", type=str, required=True)
    readme_update.add_argument("--path", type=str, required=False)

    metadata_update.add_argument("--repo", type=str, required=True)
    metadata_update.add_argument("--path", type=str, required=False)

    args = parser.parse_args()

    if args.command == "airtable-insert":
        print("Inserting metadata to AirTable")
        insert_metadata_to_airtable(args.model, args.contributor, args.api_key)

    elif args.command == "airtable-update":
        print("Updating metadata to AirTable")
        update_metadata_to_airtable(args.user, args.repo, args.branch, args.api_key)

    elif args.command == "readme-update":
        print("Updating README from AirTable")
        update_readme_from_airtable(args.repo, args.path)
    
    elif args.command == "metadata-update":
        print("Updating metadata from AirTable")
        update_metadata_from_airtable(args.repo, args.path)

    else:
        print("Invalid command")
        parser.print_help()
        exit(1)
