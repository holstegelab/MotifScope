import os
from setuptools import setup, find_packages

os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
here = os.getcwd()
main_script = os.path.join(here, 'motifscope')

packages = find_packages(where=here, include=['mscope', 'mscope.*'])

setup(
    name="MotifScope",
    version="1.0.0",
    packages=packages,
    scripts=[main_script],
    install_requires=[
        'numpy','pandas','matplotlib','biopython','multiprocess','umap-learn','scikit-learn','levenshtein','pyabpoa'],
    author="Y.Zhang",
    author_email="y.zhang1@amsterdamumc.nl",
    description="MotifScope is a tool to characterize and visualize motif composition of tandem repeats"
)
