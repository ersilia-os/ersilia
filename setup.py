from setuptools import setup, find_packages

def get_version_and_cmdclass(package_path):
    """Load version.py module without importing the whole package.
    Template code from miniver
    """
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("version", os.path.join(package_path, "_clean_version.py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.__version__, module.cmdclass

version, cmdclass = get_version_and_cmdclass("ersilia")

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    install_requires = f.read().splitlines()

# Filter dependencies based on names and a larger list of requirements
def _filter_requires(names, requires):
    filt_reqs = []
    for r in requires:
        for n in names:
            if n in r:
                filt_reqs += [r]
    return filt_reqs

# Slim requirements
slim = [
    "bentoml",
    "PyYAML",
    "pysmiles",
    "dockerfile-parse",
    "pygit2",
    "emoji"
]
slim_requires = _filter_requires(slim, install_requires)

# Web app requirements
webapp = slim + [
    "streamlit"
]
webapp_requires = _filter_requires(webapp, install_requires)

# Doc builder requirements
doc_builder = slim + [
    "sphinx"
]
doc_builder_requires = _filter_requires(doc_builder, install_requires)

# Test requirements
test = slim + [
    "pytest"
]
test_requires = _filter_requires(test, install_requires)

# Development requires
dev_requires = install_requires

# Define extras requires
extras_require = {
    "webapp": webapp_requires,
    "doc_builder": doc_builder_requires,
    "test": test_requires,
    "dev": dev_requires,
}

setup(
    name='ersilia',
    version=version,
    cmdclass=cmdclass,
    author='Ersilia Open Source Initiative',
    author_email='miquel@ersilia.io',
    url='https://github.com/ersilia-os/ersilia',
    description='Ersilia model hub for open source drug discovery and global health',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    python_requires=">=3.7",
    install_requires=slim_requires,
    extras_require=extras_require,
    packages=find_packages(exclude=("utilities")),
    entry_points={
        'console_scripts': [
            'ersilia=ersilia.cli:cli'
        ]
    },
    classifiers=(
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence"
    ),
    keywords='drug-discovery machine-learning ersilia open-science global-health model-hub',
    project_urls={
        "Landing page": "https://ersilia.io",
        "Models": "https://ersilia-hub.netlify.app",
        "Source Code": "https://github.com/ersilia-os/ersilia/",
    },
    include_package_data=True,
)
