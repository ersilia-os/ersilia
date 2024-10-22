import os
import requests
from dotenv import load_dotenv

load_dotenv()

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
AWS_ACCESS_KEY_ID = os.getenv("AWS_ACCESS")


PRIMARY_SYSTEM_PROMPT = """
You are a biomedical expert and you are asked to annotate metadata for a given computational tool, for example an AI/ML model.
You will be given the following information:
1. A structured report (summary) of a publication. This will be provided as an attachment.
2. Some raw and potentially spurious metadata copied inline as part of the user prompt.

Your task is to provide a structured metadata report for the computational tool based on the information provided.
The metadata report should include the following information, in Markdown format:

# General Information
- Title: Suggest a title for the computational tool. This should be a concise and informative title. It should be coherent with the title of the publication, and it can be inspired by the title in the raw metadata. The title should not be longer than 100 characters.
- Slug: Just take the slug from the publication summary.

# Description
Write a short summary of the computational tool. This should be a high-level overview of the tool based on the publication summary provided.
You an also get inspiration from the user prompt.
The description should have at least 200 characters. Strictly one paragraph.

# Input
- Input type:
- Input format:

# Output
- Output type:
- Output format:


# Publication and Code

Provide a brief summary of the publication. You can use the TLDR from the publication report and the publication summary.
This should be a high-level overview of the publication. Between 50 and 100 words. Only one paragraph. No new-line characters.

In addition, provide the following URLs as a list:
- Publication URL: Provide the URL of the publication. This is exactly the URL provided in the metadata from the user.
- Code URL: Provide the URL of the code repository. This is exactly the URL provided in the metadata from the user.

"""

PRIMARY_USER_PROMPT = """


"""



class ContentFetcher(object):
    def __init__(self, model_id):
        self.model_id = model_id

    def get_publication_summary(self):
        pass

    def get_metadata_from_issue(self):
        pass


class MetadataLLM(object):
    def __init__(self, model_id):
        self.model_id = model_id

    