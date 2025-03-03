
import csv
import os

import pyairtable
import requests
import yaml

from io import StringIO

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

    def write_information_0(self, data: BaseInformation, readme_path=None):
        """
        Generates the README file using the old format.
        
        """
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

    @staticmethod
    def format_link(label, url):
        if url and isinstance(url, str) and url.startswith("http"):
            return f"[{label}]({url})"
        return label

    def write_information_1(self, data: BaseInformation, readme_path=None):
        """
        Generates the README file using the new format.
        """
        
        d = data.as_dict()
        text = "# {0}\n\n".format(d.get("Title", "Model Title"))
        
        # Description and incorporation date
        text += "**Description**  \n"
        text += "{0}\n\n".format(d.get("Description", "No description provided").rstrip("\n"))
        if "Incorporation Date" in d:
            text += "**Model Incorporation Date:** {0}\n\n".format(d["Incorporation Date"])
        text += "---\n\n"

        # Identifiers
        text += "## Identifiers\n"
        text += "- **Ersilia Identifier:** {0}\n".format(d.get("Identifier", ""))
        text += "- **Slug:** {0}\n\n".format(d.get("Slug", ""))

        # Domain
        text += "## Domain\n"
        text += "- **Task:** {0}\n".format(d.get("Task", ""))
        if d.get("Subtask"):
            text += "- **Subtask:** {0}\n".format(d.get("Subtask"))
        if d.get("Biomedical Area"):
            text += "- **Biomedical Area:** {0}\n".format(d.get("Biomedical Area"))
        if d.get("Target Organism"):
            text += "- **Target Organism:** {0}\n".format(d.get("Target Organism"))
        if d.get("Tag"):
            text += "- **Tags:** {0}\n\n".format(d.get("Tag"))
        # Input
        text += "## Input\n"
        text += "- **Input Type:** {0}\n".format(d.get("Input", ""))
        text += "- **Input Shape:** {0}\n\n".format(d.get("Input Shape", ""))

        # Output
        text += "## Output\n"
        text += "- **Output Type:** {0}\n".format(d.get("Output Type", ""))
        if d.get("Output Dimension"):
            text += "- **Output Dimension:** {0}\n".format(d.get("Output Dimension"))
        if d.get("Output Consistency"):
            text += "- **Output Consistency:** {0}\n".format(d.get("Output Consistency"))
        text += "- **Interpretation:** {0}\n\n".format(d.get("Interpretation", ""))
        # Attempt to read output columns if not present in metadata
        if not d.get("Output Columns"):
            repo_name = d.get("Identifier")
            if repo_name:
                # Construct the raw URL to the CSV file
                url = f"https://raw.githubusercontent.com/ersilia-os/{repo_name}/main/model/framework/columns/run_columns.csv"
                print("Fetching output columns from:", url)
                try:
                    response = requests.get(url)
                    #print("Response status code:", response.status_code)
                    if response.status_code == 200:
                        f = StringIO(response.text)
                        reader = csv.DictReader(f)
                        output_columns = [row for i, row in enumerate(reader) if i < 10]
                        d["Output Columns"] = output_columns
                    else:
                        d["Output Columns"] = []
                except Exception as e:
                    d["Output Columns"] = []
            else:
                d["Output Columns"] = []

        # Generate the output columns table if available
        if d.get("Output Columns"):
            text += "**Output Columns (up to 10):**\n\n"
            text += "| Name | Type | Direction | Description |\n"
            text += "|------|------|-----------|-------------|\n"
            for col in d.get("Output Columns", []):
                text += "| {0} | {1} | {2} | {3} |\n".format(
                    col.get("name", ""),
                    col.get("type", ""),
                    col.get("direction", ""),
                    col.get("description", "")
                )
            text += "\n"

        # Source and Deployment
        text += "## Source and Deployment\n"
        if d.get("Source"):
            text += "- **Source:** {0}\n".format(d.get("Source"))
        if d.get("Source Type"):
            text += "- **Source Type:** {0}\n".format(d.get("Source Type"))
        if d.get("DockerHub"):
            text += "- **DockerHub:** {0}\n".format(self.format_link("DockerHub", d["DockerHub"]))
        if d.get("Data Architecture"):
            text += "- **Data Architecture:** {0}\n".format(d.get("Data Architecture"))
        if d.get("S3"):
            # In the Source and Deployment section:
            text += "- **S3:** {0}\n".format(self.format_link("S3", d.get("S3", "")))
        if d.get("Host URL"):
            text += "- **Host URL:** {0}\n".format(self.format_link("Host URL", d["Host URL"]))
        text += "\n"

        # Resource Consumption
        text += "## Resource Consumption\n"
        if d.get("Model Size"):
            text += "- **Model Size:** {0}\n".format(d.get("Model Size"))
        if d.get("Environment Size"):
            text += "- **Environment Size:** {0}\n".format(d.get("Environment Size"))
        if d.get("Image Size"):
            text += "- **Image Size:** {0}\n".format(d.get("Image Size"))
        text += "\n"
        text += "**Computational Performance:**\n"
        if d.get("Computational Performance 1"):
            text += "- **1 input:** {0}\n".format(d.get("Computational Performance 1"))
        if d.get("Computational Performance 10"):
            text += "- **10 inputs:** {0}\n".format(d.get("Computational Performance 10"))
        if d.get("Computational Performance 100"):
            text += "- **100 inputs:** {0}\n".format(d.get("Computational Performance 100"))
        text += "\n"

        text += "## Additional Information\n"
        if d.get("DO Deployment"):
            text += "- **DO Deployment:** {0}\n".format(d.get("DO Deployment"))
        if d.get("Biomodel Annotation"):
            text += "- **Biomodel Annotation:** {0}\n".format(d.get("Biomodel Annotation"))
        if d.get("Runtime"):
            text += "- **Runtime:** {0}\n".format(d.get("Runtime"))
        if d.get("Secrets"):
            text += "- **Secrets:** {0}\n".format(d.get("Secrets"))
        if d.get("Deployment"):
            text += "- **Deployment:** {0}\n".format(d.get("Deployment"))
        if d.get("Incorporation Quarter"):
            text += "- **Incorporation Quarter:** {0}\n".format(d.get("Incorporation Quarter"))
        if d.get("Incorporation Year"):
            text += "- **Incorporation Year:** {0}\n".format(d.get("Incorporation Year"))
        text += "\n"

        # References
        text += "## References\n"
        text += "- **Source Code:** {0}\n".format(self.format_link("Source Code", d.get("Source Code", "")))
        text += "- **Publication:** {0}\n".format(self.format_link("Publication", d.get("Publication", "")))
        if "Publication Type" in d:
            text += "  - **Publication Type:** {0}\n".format(d["Publication Type"])
        if "Publication Year" in d:
            text += "  - **Publication Year:** {0}\n".format(d["Publication Year"])
        text += "\n"

        # License
        text += "## License\n"
        license_value = d.get("License", "").strip()
        if license_value:
            license_text = (f"This package is licensed under a GPL-3.0 license. "
                            f"The model contained within this package is licensed under a {license_value} license.")
        else:
            license_text = ("This package is licensed under a GPL-3.0 license. "
                            "The model contained within this package is licensed under a MIT license.")
        text += license_text + "\n\n"
        text += ("Notice: Ersilia grants access to these models 'as is' provided by the original authors, "
                "please refer to the original code repository and/or publication if you use the model in your research.\n\n")

        # About Ersilia
        text += "## About Ersilia\n"
        text += "The [Ersilia Open Source Initiative](https://ersilia.io) is a non-profit organization fueling sustainable research in the Global South.\n\n"
        text += "[Help us](https://www.ersilia.io/donate) achieve our mission!"
        if readme_path is None:
            return text
        else:
            with open(readme_path, "w") as f:
                f.write(text)

    def write_information(self, data: BaseInformation, readme_path=None, version=None):
        """
        Dynamically selects the appropriate README format based on the provided version flag
        or by inspecting the metadata.
        
        :param data: A BaseInformation instance containing the model metadata.
        :param readme_path: If provided, the README will be written to this path.
        :param version: '0' for old format, '1' for new format. If None, auto-detection is attempted.
        :return: The generated README text (if readme_path is None).
        """
        d = data.as_dict()
        if version is None:
            # Auto-detect based on the presence of keys specific to the new format.
            if "Output Columns" in d or "Model Incorporation Date" in d:
                version = 1
            else:
                version = 0
        if version == 0:
            return self.write_information_0(data, readme_path)
        elif version == 1:
            return self.write_information_1(data, readme_path)
        else:
            raise ValueError("Unknown version specified for README generation.")


class ReadmeUpdater:
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
    readme_update.add_argument("--path", type=str, required=True)

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

    else:
        print("Invalid command")
        parser.print_help()
        exit(1)
