import sys
from setuptools import setup, find_packages

packages = find_packages(include=['mscope', 'mscope.*'])

setup(
    name="MotifScope",
    version="0.1.0",
    packages = packages,
    scripts = ['motifscope'],
    install_requires=['numpy','pandas','matplotlib','biopython', 'multiprocess', 'umap-learn', 'scikit-learn', 'levenshtein', 'pyabpoa'],
    author = "Y.Zhang",
    author_email = "y.zhang1@amsterdamumc.nl",
    description = "MotifScope is a tool to characterize and visualize motif composition of tandem repeats"
)

