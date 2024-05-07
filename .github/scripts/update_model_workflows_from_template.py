import sys
import os
import shutil
import subprocess

model_repo = sys.argv[1]

template_repo = "eos-template"

subprocess.Popen(
    "git clone https://github.com/ersilia-os/eos-template.git", shell=True
).wait()
subprocess.Popen(
    "git clone https://github.com/ersilia-os/{0}.git".format(model_repo), shell=True
).wait()

github_folder = os.path.join("{0}/.github/".format(model_repo))

if os.path.exists(github_folder):
    shutil.rmtree(github_folder)

shutil.copytree("eos-template/.github/", github_folder)

subprocess.Popen(
    "cd {0};git add .; git commit -m 'updated workflow files'; git push; cd ..".format(
        model_repo
    ),
    shell=True,
).wait()

subprocess.Popen("rm -rf eos-template", shell=True).wait()
subprocess.Popen("rm -rf {0}".format(model_repo), shell=True).wait()
