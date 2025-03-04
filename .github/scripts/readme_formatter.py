import csv
import requests


class ReadmeFormatter():
    def __init__(self):
        pass
    
    def write_information_0(self, d):
        """
        Generates the README file using the old format.
        
        """
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
        text += "* [GitHub]({0})\n".format(d["GitHub"])
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
        return text

    @staticmethod
    def format_link(label, url):
        if url and isinstance(url, str) and url.startswith("http"):
            return f"[{label}]({url})"
        return label

    def write_information_1(self, d):
        """
        Generates the README file using the new format.
        """
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
                    if response.status_code == 200:
                        reader = csv.DictReader(response.text.splitlines())
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
        if d.get("Docker Architecture"):
            text += "- **Docker Architecture:** {0}\n".format(d.get("DDocker Architecture"))
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

        # References
        text += "## References\n"
        text += "- **Source Code:** {0}\n".format(self.format_link("Source Code", d.get("Source Code", "")))
        text += "- **Publication:** {0}\n".format(self.format_link("Publication", d.get("Publication", "")))
        if "Publication Type" in d:
            text += "  - **Publication Type:** {0}\n".format(d["Publication Type"])
        if "Publication Year" in d:
            text += "  - **Publication Year:** {0}\n".format(d["Publication Year"])
        text += "\n"
        if d.get("Incorporation Date"):
            text += "- **Incorporation Date:** {0}\n".format(d.get("Incorporation Date"))    

        # License
        text += "## License\n"
        license_value = d.get("License", "").strip()
        if license_value:
            license_text = (f"This package is licensed under a GPL-3.0 license. "
                            f"The model contained within this package is licensed under a {license_value} license.")
        else:
            license_text = ("This package is licensed under a GPL-3.0 license.")
        text += license_text + "\n\n"
        text += ("Notice: Ersilia grants access to these models 'as is' provided by the original authors, "
                "please refer to the original code repository and/or publication if you use the model in your research.\n\n")

        # About Ersilia
        text += "## About Ersilia\n"
        text += "The [Ersilia Open Source Initiative](https://ersilia.io) is a non-profit organization fueling sustainable research in the Global South.\n\n"
        text += "[Help us](https://www.ersilia.io/donate) achieve our mission!"
        return text