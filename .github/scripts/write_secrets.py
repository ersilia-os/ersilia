import collections
import json
import os

OUTPUT = "secrets.json"

PATH = os.path.dirname(os.path.realpath(__file__))

PHRASE = "Write secrets to json file"

secrets = []
with open(os.path.join(PATH, "../workflows/write_secrets.yml"), "r") as f:
    text = f.read()
    chunks = text.split("- name: ")
    mychunk = None
    for chunk in chunks:
        if chunk[: len(PHRASE)] == PHRASE:
            mychunk = chunk
    envs = mychunk.split("env:\n")[1].split("run:")[0].split("\n")
    for env in envs:
        if "secrets." in env:
            secrets += [env.strip().split(":")[0]]

data = collections.OrderedDict()
for secret in secrets:
    data[secret] = os.environ[secret]

with open(os.path.join(PATH, "../..", OUTPUT), "w") as f:
    json.dump(data, f, indent=4)
