from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    description = fh.read()

setup(
    name="Clusttree",
    version="0.1.0",
    author="Yann Sevellec",
    author_email="y.sevellec@live.fr",
    packages=find_packages(),
    description="a python package for automated generation of population genomic phylogenetic trees from genbank files, it include automated closest reference detection, clustering along the references and generation of population level variant claling and variant annotation and generation of a tree of life like community level phylogenetic tree.",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/ysevel/Clusttree",
    license='incoming',
    python_requires='=3.9',
    install_requires=[
                      'snippy>=4.6',
                      'bakta>=11.4'],
    entry_points={'console_scripts': ["Clusttree=src.phylogeny_workflow:main"]}
)
