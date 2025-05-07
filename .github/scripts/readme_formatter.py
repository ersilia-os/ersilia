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
        text += "* Computational Performance For Ten Inputs: `{0}`\n".format(d["Computational Performance 1"])
        text += "* Computational Performance For Hundred Inputs: `{0}`\n".format(d["Computational Performance 3"])
        text += "* Computational Performance For Ten Thousand Inputs: `{0}`\n".format(d["Computational Performance 5"])
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
        text += "{0}\n\n".format(d.get("Description", "No description provided").rstrip("\n"))
        if "Incorporation Date" in d:
            text += "This model was incorporated on {0}.\n".format(d.get("Incorporation Date"))
        text += "\n"

        # Information section
        text += "## Information\n"

        # Identifiers
        text += "### Identifiers\n"
        text += "- **Ersilia Identifier:** `{0}`\n".format(d.get("Identifier", ""))
        text += "- **Slug:** `{0}`\n".format(d.get("Slug", ""))
        text += "\n"

        # Domain
        text += "### Domain\n"
        text += "- **Task:** `{0}`\n".format(d.get("Task", ""))
        if d.get("Subtask"):
            text += "- **Subtask:** `{0}`\n".format(d.get("Subtask"))
        if d.get("Biomedical Area"):
            text += "- **Biomedical Area:** {0}\n".format(", ".join(["`{0}`".format(x) for x in d.get("Biomedical Area")]))
        if d.get("Target Organism"):
            text += "- **Target Organism:** {0}\n".format(", ".join(["`{0}`".format(x) for x in d.get("Target Organism")]))
        if d.get("Tag"):
            text += "- **Tags:** {0}\n".format(", ".join(["`{0}`".format(x) for x in d.get("Tag")]))
        text += "\n"
        
        # Input
        text += "### Input\n"
        if d.get("Input"):
            text += "- **Input:** {0}\n".format(", ".join(["`{0}`".format(x) for x in d.get("Input", "")]))
        if d.get("Input Dimension"):
            text += "- **Input Dimension:** `{0}`\n".format(d.get("Input Dimension", ""))
        text += "\n"

        # Output
        text += "### Output\n"
        if d.get("Output"):
            output_type = d.get("Output Type", "")
            if output_type:
                text += "- **Output:** `{0}`\n".format(d.get("Output Type", ""))
        if d.get("Output Dimension"):
            text += "- **Output Dimension:** `{0}`\n".format(d.get("Output Dimension"))
        if d.get("Output Consistency"):
            text += "- **Output Consistency:** `{0}`\n".format(d.get("Output Consistency"))
        text += "- **Interpretation:** {0}\n\n".format(d.get("Interpretation", ""))

        repo_name = d.get("Identifier")
        if repo_name:
            # Construct the raw URL to the CSV file
            url = f"https://raw.githubusercontent.com/ersilia-os/{repo_name}/main/model/framework/columns/run_columns.csv"
            print("Fetching output columns from:", url)
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    reader = csv.DictReader(response.text.splitlines())
                    output_columns = [row for i,row in enumerate(reader)]
                    output_columns_short = output_columns[:10]
                    d["Output Columns"] = output_columns_short
                else:
                    d["Output Columns"] = []
            except Exception as e:
                d["Output Columns"] = []
        else:
            d["Output Columns"] = []

        # Generate the output columns table if available
        if d.get("Output Columns"):
            text += "Below are the **Output Columns** of the model:\n"
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
            if len(output_columns)>10:
                text += f"_10 of {len(output_columns)} columns are shown_"
            text += "\n"

        # Source and Deployment
        text += "### Source and Deployment\n"
        if d.get("Source"):
            text += "- **Source:** `{0}`\n".format(d.get("Source"))
        if d.get("Source Type"):
            text += "- **Source Type:** `{0}`\n".format(d.get("Source Type"))
        if d.get("DockerHub"):
            text += "- **DockerHub**: [{0}]({0})\n".format(d.get("DockerHub"))
        if d.get("Docker Architecture"):
            text += "- **Docker Architecture:** {0}\n".format(", ".join(["`{0}`".format(x) for x in d.get("Docker Architecture")]))
        if d.get("S3"):
            text += "- **S3 Storage**: [{0}]({0})\n".format(d.get("S3"))
        if d.get("Host URL"):
            text += "- **Host URL**: [{0}]({0})\n".format(d.get("Host URL"))
        text += "\n"

        # Resource Consumption
        text += "### Resource Consumption\n"
        if d.get("Model Size"): #assuming Model Size will be the first field filled here
            text += "- **Model Size (Mb):** `{0}`\n".format(d.get("Model Size"))
        if d.get("Environment Size"):
            text += "- **Environment Size (Mb):** `{0}`\n".format(d.get("Environment Size"))
        if d.get("Image Size"):
            text += "- **Image Size (Mb):** `{0}`\n".format(d.get("Image Size"))
        text += "\n"
        if d.get("Computational Performance 1"):
            text += "**Computational Performance (seconds):**\n"
        if d.get("Computational Performance 1"):
            text += "- 10 inputs: `{0}`\n".format(d.get("Computational Performance 1"))
        if d.get("Computational Performance 3"):
            text += "- 100 inputs: `{0}`\n".format(d.get("Computational Performance 3"))
        if d.get("Computational Performance 5"):
            text += "- 10000 inputs: `{0}`\n".format(d.get("Computational Performance 5"))
        text += "\n"

        # References
        text += "### References\n"
        text += "- **Source Code**: [{0}]({0})\n".format(d.get("Source Code"))
        text += "- **Publication**: [{0}]({0})\n".format(d.get("Publication"))
        if "Publication Type" in d:
            text += "- **Publication Type:** `{0}`\n".format(d.get("Publication Type"))
        if "Publication Year" in d:
            text += "- **Publication Year:** `{0}`\n".format(d.get("Publication Year"))
        if d.get("Contributor"):
            text += "- **Ersilia Contributor:** [{0}](https://github.com/{0})\n".format(
            d["Contributor"]
        )   
        text += "\n"

        # License
        text += "### License\n"
        license_value = d.get("License", "").strip()
        if license_value:
            license_text = (f"This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. "
                            f"The model contained within this package is licensed under a [{license_value}](LICENSE) license.")
        else:
            license_text = ("This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license.")
        text += license_text + "\n\n"
        text += ("**Notice**: Ersilia grants access to models _as is_, directly from the original authors, "
                "please refer to the original code repository and/or publication if you use the model in your research.\n")
        text += "\n\n"

        # How to use the model
        text += "## Use\n"
        text += "To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.\n"
        text += "The model can be **fetched** using the following command:\n"
        text += "```bash\n"
        text += "# fetch model from the Ersilia Model Hub\n"
        text += f"ersilia fetch {d.get('Identifier')}\n"
        text += "```\n"
        text += "Then, you can **serve**, **run** and **close** the model as follows:\n"
        text += "```bash\n"
        text += "# serve the model\n"
        text += f"ersilia serve {d.get('Identifier')}\n"
        text += "# generate an example file\n"
        text += f"ersilia example -n 3 -f my_input.csv\n"
        text += "# run the model\n"
        text += f"ersilia run -i my_input.csv -o my_output.csv\n"
        text += "# close the model\n"
        text += f"ersilia close\n"
        text += "```\n\n"
        
        # About Ersilia
        text += "## About Ersilia\n"
        text += "The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.\n"
        text += "Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.\n"
        text += "If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!\n"
        return text