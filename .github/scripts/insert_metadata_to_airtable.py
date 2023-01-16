import sys
import json
from pyairtable import Table
import requests

write_api_key = sys.argv[2]
model_id = sys.argv[1]

base_id = "appgxpCzCDNyGjWc8"
table_name = "Models"
METADATA_JSON_FILE = "metadata.json"
org = "ersilia-os"
branch = "main"
json_path = "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
    org, model_id, branch, METADATA_JSON_FILE
)
github = "https://github.com/ersilia-os/{0}".format(model_id)

r = requests.get(json_path)
if r.status_code == 200:
    text = r.content
    data = json.loads(text)

data["GitHub"] = github
data["Status"] = "In Progress"

table = Table(write_api_key, base_id, table_name)
table.create(data)
