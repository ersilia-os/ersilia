import os
import sys
import requests
import shutil

BUILD_VERSIONS = ["ersiliapack", "legacy-bentoml", "multistage-condapack"]
ENV_TYPES = ["conda", "pip"]
model_id = sys.argv[1]
build_version = (
    sys.argv[2] if len(sys.argv) > 2 else "multistage-condapack"
)  # This is the most stable right now.
env_type = (
    "." + sys.argv[3] if (len(sys.argv) > 3 and sys.argv[2] == "ersiliapack") else ""
)

if os.path.exists("Dockerfile"):
    shutil.move("Dockerfile", "Dockerfile_legacy")


def download_file(url, filename):
    r = requests.get(url, allow_redirects=True)
    open(filename, "wb").write(r.content)


url = f"https://raw.githubusercontent.com/ersilia-os/ersilia/master/dockerfiles/dockerize-{build_version}/model/Dockerfile{env_type}"
print("URL: ", url)
filename = "Dockerfile"
download_file(url, filename)

with open("Dockerfile", "r") as f:
    text = f.read()

text = text.replace("eos_identifier", model_id)

with open("Dockerfile", "w") as f:
    f.write(text)
