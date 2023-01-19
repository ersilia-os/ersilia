import sys
from ersilia.publish.metadata import ReadmeUpdater

repo_name = sys.argv[1]

rm = ReadmeUpdater(model_id=repo_name)
rm.update()
