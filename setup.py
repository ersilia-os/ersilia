from setuptools import setup, find_packages 

with open("README.md", "r") as fh:
    long_description = fh.read()
  
setup( 
        name ='modelhub-ai', 
        version ='0.0.6', 
        author ='Ahmed Hosny', 
        author_email ='info@modelhub.ai', 
        url ='https://github.com/modelhub-ai/modelhub', 
        description ='', 
        long_description = long_description, 
        long_description_content_type ="text/markdown", 
        license ='MIT', 
        packages = find_packages(exclude=("utilities", "examples")), 
        entry_points ={ 
            'console_scripts': [ 
                'modelhub = modelhub.start:main',
                'modelhub-ai = modelhub.start:main',
                'modelhub-run = modelhub.start:main',
                'modelhub-list = modelhub.start:list_online_models'
            ] 
        }, 
        classifiers =( 
            "Programming Language :: Python :: 3", 
            "License :: OSI Approved :: MIT License", 
            "Operating System :: OS Independent", 
        ), 
        keywords ='ai modelhub harvard dana-farber brigham deep learning open science',
        zip_safe = False
) 