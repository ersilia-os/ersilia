import os
import argparse
import tempfile
import PyPDF2
import openai
import requests
import boto3
import shutil
from dotenv import load_dotenv

load_dotenv()

# Authenticate into OpenAI
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
openai.api_key = OPENAI_API_KEY
MODEL_NAME = "gpt-4o"

# Authenticate into AWS
aws_access_key_id = os.getenv("AWS_ACCESS_KEY")
aws_secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")

if aws_access_key_id and aws_secret_access_key:
    boto3.setup_default_session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
    )
else:
    raise Exception("AWS credentials are not set in the environment variables")


PRIMARY_SYSTEM_PROMPT = """
Your are a biomedical expert. You have to make a structured report of a scientific publication. The publication will be provided as a PDF file.
Your report needs to be concise and informative. Strictly follow this structure, in Markdown format:

# Publication details
- Title: Just copy the title of the publication
- Authors: List the authors separated by commas. For example, G. Turon, M. Duran-Frigola and D. Arora.
- Journal: Name of the journal
- Year: Year of publication
- Suggested slug: A short version of the title, with hyphens instead of spaces. For example, "deep-learning-for-malaria". The slug should be all lowercase. It cannot be longer than 50 characters and should not contain more than 5 hyphens. If there is a name for the method, for example, chemprop, try to use it in the slug. Do not add dates or names of authors. If possible, do not include words such as ml (for machine learning) or ai (for artificial intelligence).
- Suggested computational title: Suggest a title that is focused on the computational methods used in the publication. For example, "Broad-spectrum antibiotics activity prediction" or "MAIP, antimalarial activity prediction based on multiple industry datasets".

# TLDR
Write a short summary of the publication in one or two sentences. This should be a high-level overview of the publication.
Do not use new-line characters. The TLDR should have between 100 and 200 characters.

# Summary
Write a summary of the publication. Feel free to use the information on the Abstract, if available.
The summary should be between 100 and 200 words. Only one paragraph. No new-line characters are allowed. No special characters, boldface, links or references are allowed.
Use a concise style.

# Relevance to biomedical research
Briefly discuss why the publication is relevant to biomedicine or drug discovery.
If the publication is related to a particular disease or pathogen, make sure to mention it.
If the publication is related to a particular stage of the drug discovery pipeline, make sure to mention it.
This should be between 50 and 100 words. Only one paragraph. No new-line characters.
Use a concise style.

# Computational methods
Write a summary of the computational methods used in this publication.
If AI/ML methods are used, focus on those. Mention the main techniques and methods.
Try to explain what the input and output of the methods are. In addition, try to explain how to interpret the output and what range of values are expected and relevant (if applicable).
If training data was used, mention the size and source of the data.
If accuracy was reported, make sure to mention it.
The summary should be between 100 and 200 words. Only one paragraph. No new-line characters are allowed. No special characters, boldface, links or references are allowed.
Use a concise style.

# Biomedical keywords
Suggest 3 to 5 keywords that are relevant to the publication. These keywords should be related to the biomedical or drug discovery aspects of the publication.
If the publication is not related to biomedicine, do not suggest any keywords.
- Keyword 1
- Keyword 2
- Keyword 3

# Computational keywords
Suggest 3 to 5 keywords that are relevant to the publication. These keywords should be related to the computational aspects of the publication, especially the AI/ML aspects, if applicable.
Do not use keywords that are too generic, such as "machine learning", "deep learning" or "descriptors".
- Keyword 1
- Keyword 2
- Keyword 3

# Strenghts
Discuss the strengths of the publication, especially from the perspective of the computational methods and the training dataset, if applicable.
Why are the results of the publication important or relevant? What are the main contributions of the publication?
This should be between 50 and 100 words. Only one paragraph. No new-line characters.

# Limitations
Discuss the limitations of the publication, especially from the perspective of the computational methods and the training dataset, if applicable.
What could be improved in the publication? What are the main weaknesses? Are the computational methods novel? Are the results reliable? Are the conclusions valid?
Is the dataset large enough? Is the data of high quality?
This should be between 50 and 100 words. Only one paragraph. No new-line characters.

# Overall relevance
Try to assess the relevance of the publication in the context of the current knowledge in the field.
In your assessment, consider the novelty of the methods, the quality of the results, and the potential impact of the publication.
The date of publication is also important.
The size of the dataset and the quality of the data are also important factors to consider.
If prior art on the topic or similar articles exist in the literature, this should penalize the relevance. Novelty is important.
Also consider the performance of the computational methods, and compare it with the performance of other methods.
The impact factor of the journal is an important factor for relevance. Higher relevance should be given to higher impact journals.
Do not be over-emphatic. Try to be explicit about the high, medium or low relevance of the publication. Not all publications are highly relevant, so do a fair assessment.
This should be between 50 and 100 words. Only one paragraph. No new-line characters. Be concise.

---

Below are some style guidelines for the report:
- Always use the third person and do not use personal pronouns.
- Do not start paragraphs with sentences such as "this study", "this publication", "this report" or "the authors". Directly explain the content. For example, do not say "This study develops a method for...". Instead, say "A method was developed for..." or similar.
- Do not include references or links.
- Do not include any special characters, boldface, or italics.
- Each section should be a separate paragraph (one and only one paragraph per section) or a bullet point list, as applicable.
- Do not include any information about the number of pages of the publication.
- Do not mention the figures or tables in the publication.
- Do not mention other references in the publication.
- Do not mention funding sources or acknowledgements.
- Do not begin your answer with a sentence like: "here is a structured report of the publication". Start with the report itself.
- Do not end your answer with a conclusion or a summary. The report should be self-contained.
"""

PRIMARY_USER_PROMPT = """
Make a structured report of the scientific publication from the PDF file. The following is the extracted text from the PDF: 
"""

POSTPROCESS_SYSTEM_PROMPT = """
You have to standardize the format of a structured report that may contain formatting errors or might require slight modifications.
You will be provided with a text file with the report. You need to make sure that the report follows the correct format. Do not make drastic changes.
The format is as follows:

# Publication details
- Title: Title of the publication
- Authors: List of authors. For example, G. Turon, M. Duran-Frigola and D. Arora. Always abbreviate the given names.
- Journal: Name of the journal. For example, Nature Biotechnology.
- Year: Year of publication. For example, 2024.
- Suggested slug: lowercase, maximum 5 hyphens and maxiumum 50 characters. No dates or author names. If possible, no generic abbreviations such as "ml" or "ai".

# TLDR
One paragraph only. Between 100 and 200 characters.

# Summary
One paragraph only. Between 100 and 200 words.

# Relevance to biomedical research
One paragraph only. Between 50 and 100 words.

# Computational methods
One paragraph only. Between 100 and 200 words.

# Biomedical keywords
Keywords as a list and without repeats or redundancy. Minim 3 and maximum 5 keywords. For example:
- Malaria
- Plasmodium falciparum
- Asexual blood stage

# Computational keywords
Keywords as a list and without repeats or redundancy. Minim 3 and maximum 5 keywords. For example:
- Physicochemical descriptors
- Feature selection
- Support vector machines

# Strenghts
One paragraph only. Between 50 and 100 words.

# Limitations
One paragraph only. Between 50 and 100 words.

# Overall relevance
One paragraph only. Between 50 and 100 words.

---

Below are a few general rules:
- Rephrase if necessary if the content is not clear or concise, but do not make big changes. Just make the small necessary adjustments to follow the guidelines.
- Do not include any special characters, boldface, or italics.
- Do not include any references or links.
- The report should be self-contained and not require additional information.
- It is a Markdown file, so make sure to follow the Markdown syntax. All sections are title sections (#) and the lists are bullet points (-).
- Avoid using double spaces or double new lines. Use only one space or one new line.
- Correct any inconsistency. For example, if the model is deemed to be highly relevant in a section, and of medium relevance elsewhere, harmonize this. Be conservative.

Make sure that the report follows the correct format. Do not make drastic changes.
"""

POSTPROCESS_USER_PROMPT = """
Reformat, if necessary, the following publication report:
"""


class PublicationPDFDownloader(object):
    def __init__(self, model_id=None, issue_number=None, url=None, file_path=None):
        if (
            model_id is None
            and issue_number is None
            and url is None
            and file_path is None
        ):
            raise Exception("At least one argument is necessary")
        if model_id is not None:
            if issue_number is not None or url is not None or file_path is not None:
                raise Exception("Only one argument is accepted")
        if issue_number is not None:
            if model_id is not None or url is not None or file_path is not None:
                raise Exception("Only one argument is accepted")
        if url is not None:
            if (
                model_id is not None
                or issue_number is not None
                or file_path is not None
            ):
                raise Exception("Only one argument is accepted")
        if file_path is not None:
            if model_id is not None or issue_number is not None or url is not None:
                raise Exception("Only one argument is accepted")
        self.model_id = model_id
        self.issue_number = issue_number
        self.url = url
        self.file_path = file_path
        self.tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        self.output_file = os.path.join(self.tmp_folder, "publication.pdf")

    def download_by_model_id(self):
        s3 = boto3.client("s3")
        bucket_name = "model-publications"
        key = f"{self.model_id}_publication.pdf"
        try:
            with open(self.output_file, "wb") as f:
                s3.download_fileobj(bucket_name, key, f)
            return self.output_file
        except Exception as e:
            raise Exception(f"Failed to download file from S3: {str(e)}")

    def download_by_issue_number(self):
        # TODO download publication PDF from a TBD repository where publications are stored by issue number
        pass

    def download_from_url(self):
        response = requests.get(self.url)
        if response.status_code == 200:
            with open(self.output_file, "wb") as f:
                f.write(response.content)
            return self.output_file
        else:
            raise Exception(
                f"Failed to download file from URL: {self.url}, status code: {response.status_code}"
            )

    def download_from_file_path(self):
        shutil.copy(self.file_path, self.output_file)
        return self.output_file

    def download(self):
        if self.model_id is not None:
            return self.download_by_model_id()
        if self.issue_number is not None:
            return self.download_by_issue_number()
        if self.url is not None:
            return self.download_from_url()
        if self.file_path is not None:
            return self.download_from_file_path()
        return None


class PublicationSummarizer(object):
    def __init__(self, publication_pdf, output_markdown):
        self.tmp_dir = tempfile.mkdtemp(prefix="ersilia-")
        self.model_name = MODEL_NAME
        self.publication_pdf = publication_pdf
        self.output_markdown = output_markdown

    def extract_text_from_pdf(self, pdf_path):
        with open(pdf_path, "rb") as file:
            reader = PyPDF2.PdfReader(file)
            text = ""
            for page_num in range(len(reader.pages)):
                page = reader.pages[page_num]
                text += page.extract_text()
        return text

    def make_primary_request(self, extracted_text):
        system_prompt = PRIMARY_SYSTEM_PROMPT.strip()
        user_prompt = PRIMARY_USER_PROMPT.strip() + "\n" + extracted_text
        response = openai.chat.completions.create(
            model=self.model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.7,
        )
        return response.choices[0].message.content

    def postprocess_response(self, text_response, markdown_file):
        system_prompt = POSTPROCESS_SYSTEM_PROMPT.strip()
        user_prompt = POSTPROCESS_USER_PROMPT.strip() + "\n" + text_response
        response = openai.chat.completions.create(
            model=self.model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.7,
        )
        text_response = response.choices[0].message.content
        with open(markdown_file, "w") as f:
            f.write(text_response)

    def run(self):
        text = self.extract_text_from_pdf(self.publication_pdf)
        text_response = self.make_primary_request(text)
        self.postprocess_response(text_response, self.output_markdown)


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
        "-u", "--url", type=str, default=None, help="URL of the downloadable file"
    )
    parser.add_argument(
        "-f",
        "--file_path",
        type=str,
        default=None,
        help="File path of the publication PDF",
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
    url = args.url
    file_path = args.file_path
    output_markdown = args.output_markdown
    ppd = PublicationPDFDownloader(model_id, issue_number, url, file_path)
    publication_pdf = ppd.download()
    ps = PublicationSummarizer(publication_pdf, output_markdown)
    ps.run()
