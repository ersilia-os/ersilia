import os
import argparse
import openai
import requests
import boto3
import tempfile
import shutil
from dotenv import load_dotenv

load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
openai.api_key = OPENAI_API_KEY
MODEL_NAME = "gpt-4o"


def read_accepted_labels(category):
    category = category.lower().replace(" ", "_")
    base_url = "https://raw.githubusercontent.com/ersilia-os/ersilia/refs/heads/master/ersilia/hub/content/metadata/"
    url = os.path.join(base_url, category + ".txt")
    response = requests.get(url)
    if response.status_code == 200:
        values = response.text.splitlines()
        return [x.strip() for x in values if x]
    else:
        raise Exception(f"Failed to download labels from {url}")


accepted = {}
for category in [
    "input",
    "input_shape",
    "license",
    "mode",
    "output",
    "output_shape",
    "output_type",
    "tag",
    "task",
]:
    accepted[category] = ", ".join(read_accepted_labels(category))


PRIMARY_SYSTEM_PROMPT = """
You are a biomedical expert and you are asked to annotate metadata for a given computational tool, for example an AI/ML model.
You will be given the following information:
1. A structured report (summary) of a publication. This will be labelled by the user as PUBLICATION REPORT.
2. Some raw and potentially spurious metadata in JSON, YAML or Markdown format. This will be labeled by the user as RAW METADATA.

Your task is to provide a structured metadata report for the computational tool based on the information provided.
The metadata report should include the following information, in Markdown format:

# General Information
- Title: Suggest a title for the computational tool. This should be a concise and informative title. It should be coherent with the title of the publication, and it can be inspired by the title in the raw metadata. The title should not be longer than 100 characters.
- Slug: Just take the slug from the publication summary.

# Description
Write a short summary of the computational tool. This should be a high-level overview of the tool based on the publication summary provided.
You an also get inspiration from the user prompt.
The description should have at least 200 characters. Strictly one paragraph.

# Tags
Select 1-5 tags from this list: {7}
Give the tags as a list:
- Tag 1
- Tag 2
- etc.

# Input
- Input: One of the following, depending on the type of input accepted for the model: {0}
- Input shape: One of the following, depending on the characteristics of the minimum model input. Typically, models accept a single input. Models that require multiple inputs to do predictions, can be lists, etc.: {1}

# Output
- Output: One (ideally) or multiple of the following (comma-separated), depending on the type of output produced by the model: {2}
- Output shape: One of the following, depending on the dimensionality of the produced output: {3}
- Output type: One of the following: {4}

# Mode and Task
- Mode: One of the following, depending on whether the model is already pre-trained by authors, retrained by our team, or build from the data with our own tools: {5}
- Task: AI, ML or data science task. One of the following: {6}

# Interpretation
Provide a oneliner explaining how to interpret the output of the model. Is it a probability? Is it a regression value for a particular experimental assay? etc. No more than 150 characters.

# Publication and Code
Provide a brief summary of the publication. You can use the TLDR from the publication report and the publication summary.
This should be a high-level overview of the publication. Between 50 and 100 words. Only one paragraph. No new-line characters.

In addition, provide the following URLs as a list:
- Publication URL: Provide the URL of the publication. This is exactly the URL provided in the metadata from the user.
- Code URL: Provide the URL of the code repository. This is exactly the URL provided in the metadata from the user.

# License
- Simply extract the license from the metadata.

--

Below are some general guidelines:
- Produce a Markdown file strictly following the headers and formatting specified above.
- Never do multiple paragraphs.
- Be concise.
- Do not include any special characters, boldface, or italics.
""".format(
    accepted["input"],
    accepted["input_shape"],
    accepted["output"],
    accepted["output_shape"],
    accepted["output_type"],
    accepted["mode"],
    accepted["task"],
    accepted["tag"],
)


PRIMARY_USER_PROMPT = """
Write the metadata of the following computational tool.

RAW METADATA
<RAW_METADATA>

PUBLICATION REPORT
<PUBLICATION_REPORT>
"""

POSTPROCESS_SYSTEM_PROMPT = """
Your task is to make sure the format of a Markdown file is correct. You can slightly modify the content, but only if you find
inconsistencies in the text, repetitions, or incoherences. Avoid boldface, italics and special characters.
Strictly follow the format below:
# General information
- Title: The title. No longer than 100 characters.
- Slug: Lower case, can use hyphenation.

# Description
At least 200 characters, ideally more. Strictly one paragraph.

# Tags
- Tag 1
- Tag 2
- etc. 3 to 5 tags.

# Input
- Input: As specified in the source file. Do not modify.
- Input shape: As specified in the source file. Do not modify.

# Output
- Output: As specified in the source file. Do not modify.
- Output shape: As specified in the source file. Do not modify.
- Output type: As specified in the source file. Do not modify.

# Mode and Task
- Mode: As specified in the source file. Do not modify.
- Task: As specified in the source file. Do not modify.

# Interpretation
A one-liner. No more than 150 characters. Strictly one paragraph.

# Publication and Code
A summary as specified in the source file, between 50 and 100 words. Strictly one paragraph.

- Publication: A URL with the link []()
- Code: A URL with the link []()

# License
- As specified in the source file. Do not modify.
"""

POSTPROCESS_USER_PROMPT = """
Process the following text:
"""


class MetadataDownloader(object):
    def __init__(
        self,
        model_id=None,
        issue_number=None,
    ):
        if model_id is None and issue_number is None:
            raise Exception("At least one argument is necessary")
        if model_id is not None and issue_number is not None:
            model_id = None
        self.model_id = model_id
        self.issue_number = issue_number
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.output_file = os.path.join(self.tmp_folder, "metadata.txt")

    def download_by_model_id(self):
        base_url = (
            "https://raw.githubusercontent.com/ersilia-os/{0}/refs/heads/main/".format(
                self.model_id
            )
        )
        metadata_yml = os.path.join(base_url, "metadata.yml")
        metadata_json = os.path.join(base_url, "metadata.json")
        try:
            url = metadata_yml
            response = requests.get(url)
            if response.status_code == 200:
                with open(self.output_file, "wb") as f:
                    f.write(response.content)
                return self.output_file
        except:
            pass
        try:
            url = metadata_json
            response = requests.get(url)
            if response.status_code == 200:
                with open(self.output_file, "wb") as f:
                    f.write(response.content)
                return self.output_file
        except:
            raise Exception("Metadata not found in model repository")

    def download_by_issue_number(self):
        # TODO extract metadata from the issue itself
        pass

    def download(self):
        if self.model_id is not None:
            return self.download_by_model_id()
        if self.issue_number is not None:
            return self.download_by_issue_number()
        return None


class PublicationSummaryDownloader(object):
    def __init__(self, model_id=None, file_path=None):
        if model_id is None and file_path is None:
            raise Exception("At least one argument is necessary")
        if model_id is not None and file_path is not None:
            model_id = None
        self.model_id = model_id
        self.file_path = file_path
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.output_file = os.path.join(self.tmp_folder, "input_file.txt")

    def download_by_model_id(self):
        s3 = boto3.client("s3")
        bucket_name = "publication-summaries"  # TODO check path in S3
        key = f"{self.model_id}_summary.md"
        try:
            with open(self.output_file, "wb") as f:
                s3.download_fileobj(bucket_name, key, f)
            return self.output_file
        except Exception as e:
            raise Exception(f"Failed to download file from S3: {str(e)}")

    def download_from_file_path(self):
        shutil.copy(self.file_path, self.output_file)
        return self.output_file

    def download(self):
        if self.model_id is not None:
            return self.download_by_model_id()
        if self.file_path is not None:
            return self.download_from_file_path()
        return None


class MetadataSuggestor(object):
    def __init__(self, metadata_txt, publication_markdown, output_markdown):
        self.model_name = MODEL_NAME
        self.metadata_txt = metadata_txt
        self.publication_markdown = publication_markdown
        self.output_markdown = output_markdown

    def make_primary_request(self):
        with open(self.metadata_txt, "r") as f:
            metadata_txt = f.read()
        with open(self.publication_markdown, "r") as f:
            publication_markdown = f.read()
        system_prompt = PRIMARY_SYSTEM_PROMPT
        user_prompt = PRIMARY_USER_PROMPT.strip()
        user_prompt = user_prompt.replace("<RAW_METADATA>", metadata_txt)
        user_prompt = user_prompt.replace("<PUBLICATION_REPORT>", publication_markdown)
        response = openai.chat.completions.create(
            model=self.model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.7,
        )
        return response.choices[0].message.content

    def postprocess_markdown(self, text, output_md):
        system_prompt = POSTPROCESS_SYSTEM_PROMPT.strip()
        user_prompt = POSTPROCESS_USER_PROMPT.strip() + "\n" + text
        response = openai.chat.completions.create(
            model=self.model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.7,
        )
        text = response.choices[0].message.content
        with open(output_md, "w") as f:
            f.write(text)

    def run(self):
        text = self.make_primary_request()
        self.postprocess_markdown(text, self.output_markdown)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process arguments.")
    parser.add_argument(
        "-m", "--model_id", type=str, default=None, help="Ersilia Model Hub identifier"
    )
    parser.add_argument(
        "-i",
        "--issue_number",
        type=int,
        default=None,
        help="GitHub issue number in the ersilia-os/ersilia repository",
    )
    parser.add_argument(
        "-s", "--summary_path", type=str, default=None, help="File path"
    )
    parser.add_argument(
        "-o",
        "--output_markdown",
        type=str,
        default=None,
        required=True,
        help="Output file in Markdown format",
    )
    args = parser.parse_args()
    model_id = args.model_id
    issue_number = args.issue_number
    summary_path = args.summary_path
    output_markdown = args.output_markdown
    md = MetadataDownloader(model_id, issue_number)
    metadata_txt = md.download()
    pd = PublicationSummaryDownloader(model_id, summary_path)
    publication_markdown = pd.download()
    ms = MetadataSuggestor(metadata_txt, publication_markdown, output_markdown)
    ms.run()
