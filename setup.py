from setuptools import setup, find_packages


def get_version(package_path):
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("version", os.path.join(package_path, "_version.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    version = module.get_version_from_static()
    return version


version = get_version("ersilia")


with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

# Slim requirements
slim = [
    "inputimeout",
    "emoji",
    "validators",
    "h5py",
    "loguru",
    "pyairtable<2",
    "PyYAML",
    "dockerfile-parse",
    "tqdm",
    "click",
    "docker",
]
slim_requires = slim

# Lake requirements
lake_requires = slim_requires + ["isaura==0.1"]

# Doc builder requirements
doc_builder_requires = slim + ["sphinx", "jinja2"]

# Test requirements
test_requires = slim + ["pytest", "fuzzywuzzy"]

# Define extras requires
extras_require = {
    "lake": lake_requires,
    "docs": doc_builder_requires,
    "test": test_requires,
}

setup(
    name="ersilia",
    version=version,
    author="Ersilia Open Source Initiative",
    author_email="hello@ersilia.io",
    url="https://github.com/ersilia-os/ersilia",
    description="A hub of AI/ML models for open source drug discovery and global health",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPLv3",
    python_requires=">=3.7",
    install_requires=slim_requires,
    extras_require=extras_require,
    packages=find_packages(),
    entry_points={"console_scripts": ["ersilia=ersilia.cli:cli"]},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    keywords="drug-discovery machine-learning ersilia open-science global-health model-hub infectious-diseases",
    project_urls={
        "Landing page": "https://ersilia.io",
        "Models": "https://ersilia.io/model-hub",
        "Source Code": "https://github.com/ersilia-os/ersilia/",
    },
    package_data={"ersilia": ["hub/content/metadata/*.txt"]},
    include_package_data=True,
)


# Install bentoml if necessary
def check_bentoml(package_path):
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location(
        "bentoml_requirement",
        os.path.join(package_path, "setup", "requirements", "bentoml.py"),
    )
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    req = module.BentoMLRequirement()
    if not req.is_bentoml_ersilia_version():
        req.install()
    return version


check_bentoml("ersilia")
