import os

BASE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "base")

def generate_conda_dockerfile():

    with open(os.path.join(BASE_PATH, "conda_versions.txt"), "r") as f:
        conda_versions = f.readlines()
    
    with open(os.path.join(BASE_PATH, "Dockerfile.conda"), "r") as f:
        dockerfile = f.readlines()
    
    for version in conda_versions:
        version = version.strip()
        with open(os.path.join(BASE_PATH, "Dockerfile.conda" + version), "w") as f:
            for line in dockerfile:
                f.write(line.replace("version", version))

def generate_pip_dockerfile():
    with open(os.path.join(BASE_PATH, "python_versions.txt"), "r") as f:
        pip_versions = f.readlines()

    with open(os.path.join(BASE_PATH, "Dockerfile"), "r") as f:
        dockerfile = f.readlines()
    
    for version in pip_versions:
        version = version.strip()
        with open(os.path.join(BASE_PATH, "Dockerfile" + version), "w") as f:
            for line in dockerfile:
                f.write(line.replace("version", version))

if __name__ == "__main__":
    generate_conda_dockerfile()
    generate_pip_dockerfile()