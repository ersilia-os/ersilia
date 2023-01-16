import sys
from ersilia.hub.content.card import RepoMetadataFile, AirtableMetadata

user_name = sys.argv[1]
repo_name = sys.argv[2]
branch = sys.argv[3]
airtable_api_key = sys.argv[4]

rm = RepoMetadataFile(model_id=repo_name, config_json=None)
data = rm.read_information(org=user_name, branch=branch)
print(data.as_dict())

am = AirtableMetadata(model_id=repo_name, config_json=None)
am.set_write_api_key(airtable_api_key)
am.write_information(data)
