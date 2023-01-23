import sys
from ersilia.publish.metadata import ReadmeUpdater

repo_name = sys.argv[1]
repo_path = sys.argv[2]

rm = ReadmeUpdater(model_id=repo_name, repo_path=repo_path)
rm.update()
