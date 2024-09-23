import sys
import json
import yaml
from pyairtable import Table
import requests

write_api_key = sys.argv[3]
model_id = sys.argv[1]
contributor_name = sys.argv[2]

base_id = "appgxpCzCDNyGjWc8"
table_name = "Models"
METADATA_YAML_FILE = "metadata.yml"
org = "ersilia-os"
branch = "main"
yaml_path = "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
    org, model_id, branch, METADATA_YAML_FILE
)
github = "https://github.com/ersilia-os/{0}".format(model_id)

r = requests.get(yaml_path)
if r.status_code == 200:
    text = r.content
    data = yaml.safe_load(text)

airtable_data = {}
airtable_data["Identifier"] = model_id
airtable_data["Slug"] = data["Slug"]
airtable_data["Title"] = data["Title"]
airtable_data["GitHub"] = github
airtable_data["Contributor"] = contributor_name
airtable_data["Status"] = "In progress"
if data["Publication"] != "":
    airtable_data["Publication"] = data["Publication"]
if data["Source Code"] != "":
    airtable_data["Source Code"] = data["Source Code"]

table = Table(write_api_key, base_id, table_name)
table.create(airtable_data)
