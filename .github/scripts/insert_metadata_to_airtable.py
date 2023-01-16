import sys
import json
from pyairtable import Table

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

with open(json_path, "r") as f:
    data = json.load(f)
data["GitHub"] = github
data["Status"] = "In Progress"

table = Table(write_api_key, base_id, table_name)
table.create(data)
