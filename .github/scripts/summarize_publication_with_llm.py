import os
import sys
import tempfile
import requests
from dotenv import load_dotenv

load_dotenv()

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
AWS_ACCESS_KEY_ID = os.getenv("AWS_ACCESS")

model_id = sys.argv[1]

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
Make a structured report of the scientific publication provided as a PDF file.
"""

POSTPROCESS_SYSTEM_PROMPT = """
You have to standardize the format of a structured report that may contain formatting errors or might require slight modifications.
You will be provided with a text file with the report. You need to make sure that the report follows the correct format. Do not make drastic changes.
The format is as follows:

# Publication details:
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
- Rephrase if necessary if the content is not clear or concise, but try not to make drastic changes. Just make the small necessary adjustments to follow the guidelines.
- Do not include any special characters, boldface, or italics.
- Do not include any references or links.
- The report should be self-contained and not require additional information.
- It is a Markdown file, so make sure to follow the Markdown syntax. All sections are title sections (#) and the lists are bullet points (-).
- Avoid using double spaces or double new lines. Use only one space or one new line.

Make sure that the report follows the correct format. Do not make drastic changes.

"""

POSTPROCESS_SYSTEM_PROMPT = """
Reformat, if necessary, the publication report provided as a text file.
"""

class PublicationSummarizer(object):
    def __init__(self, model_id):
        self.model_id = model_id
        self.tmp_dir = tempfile.mkdtemp()

    def download_publication_from_s3(self):
        pass

    def run(self):
        self.download_publication_from_s3()