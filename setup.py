from setuptools import setup, find_packages 
import versioneer

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

install_requires = [

]

extras_require = install_requires

setup(
    name='ersilia',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='Ersilia Open Source Initiative',
    author_email='hello@ersilia.io',
    url='https://github.com/ersilia-os/ersilia',
    description='Ersilia model hub for open source drug discovery',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    python_requires=">=3.8.5",
    install_requires=install_requires,
    extras_require=extras_require,
    packages=find_packages(exclude=("utilities", "examples")),
    entry_points={
        'console_scripts': [
            'ersilia=ersilia.cli:cli'
        ]
    },
    classifiers=(
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence"
    ),
    keywords='drug-discovery machine-learning ersilia open-science global-health model-hub',
    project_urls={
        "Landing page": "http://ersilia.io",
        "Models": "https://github.com/ersilia-os/eos-models/",
        "Source Code": "https://github.com/ersilia-os/ersilia/",
    },
    include_package_data=True,
)
