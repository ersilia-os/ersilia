from setuptools import setup, find_packages 

with open("README.md", "r") as fh:
    long_description = fh.read()
  
setup( 
        name='ersilia',
        version='0.0.1',
        author='Miquel Duran-Frigola',
        author_email='miquel@ersilia.io',
        url='https://github.com/ersilia-os/ersilia',
        description='',
        long_description=long_description,
        long_description_content_type="text/markdown",
        license='MIT',
        packages=find_packages(exclude=("utilities", "examples")),
        entry_points={
            'console_scripts': [ 
                'ersilia = modelhub.start:main',
                'ersilia-run = modelhub.start:main',
                'ersilia-list = modelhub.start:list_online_models'
            ] 
        }, 
        classifiers=(
            "Programming Language :: Python :: 3", 
            "License :: OSI Approved :: MIT License", 
            "Operating System :: OS Independent", 
        ), 
        keywords='drug discovery ai ersilia open science',
        zip_safe=False
) 