from setuptools import setup, find_packages


def get_version_and_cmdclass(package_path):
    """Load version.py module without importing the whole package.
    Template code from miniver
    """
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location(
        "version", os.path.join(package_path, "_clean_version.py")
    )
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass


version, cmdclass = get_version_and_cmdclass("ersilia")

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

# Slim requirements
slim = [
    "bentoml @ git+https://github.com/ersilia-os/bentoml-ersilia.git",
    "inputimeout",
    "emoji",
    "validators",
    "h5py",
    "loguru",
    "pyairtable",
    "PyYAML",
    "dockerfile-parse",
    "tqdm",
]
slim_requires = slim

# Lake requirements
lake_requires = slim_requires + ["isaura==0.1"]

# Doc builder requirements
doc_builder_requires = slim + ["sphinx", "jinja2"]

# Test requirements
test_requires = slim + ["pytest"]

# Define extras requires
extras_require = {
    "lake": lake_requires,
    "doc_builder": doc_builder_requires,
    "test": test_requires,
}

print(extras_require)

setup(
    name="ersilia",
    version=version,
    cmdclass=cmdclass,
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
