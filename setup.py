import os
from setuptools import setup

def get_version():
    v = "0.0.0"
    return v


setup(
        name = "hivevo",
        version = get_version(),
        author = "Fabio Zanini and Richard Neher",
        author_email = "richard.neher@unibas.ch",
        description = ("HIVevo access"),
        long_description = "This is a simple collection of classes that provide a python interface to precomputed data from the HIVEVO project",
        long_description_content_type="text/markdown",
        license = "MIT",
        keywords = "",
        url = "https://github.com/neherlab/HIVEVO_access",
        packages=['hivevo'],
        install_requires = [
            'biopython>=1.67',
            'numpy>=1.10.4',
            'pandas>=0.17.1',
            'scipy>=0.16.1'
        ],
        extras_require = {
            ':python_version >= "3.6"':['matplotlib>=2.0'],
        },
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8"
            ]
        )
