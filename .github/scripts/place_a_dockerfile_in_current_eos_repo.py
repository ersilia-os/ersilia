import os
import sys
import requests
import shutil

model_id = sys.argv[1]

if os.path.exists("Dockerfile"):
    shutil.move("Dockerfile", "Dockerfile_legacy")


def download_file(url, filename):
    r = requests.get(url, allow_redirects=True)
    open(filename, "wb").write(r.content)


url = "https://raw.githubusercontent.com/ersilia-os/ersilia/master/dockerfiles/model-deploy/model/Dockerfile"
filename = "Dockerfile"
download_file(url, filename)

with open("Dockerfile", "r") as f:
    text = f.read()

text = text.replace("eos_identifier", model_id)

with open("Dockerfile", "w") as f:
    f.write(text)
