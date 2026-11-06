from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    description = fh.read()

setup(
    name="Clusttree",
    version="1.0.0",
    author="Yann Sevellec",
    author_email="y.sevellec@live.fr",
    packages=find_packages(),
    description="a python package for automated generation of population genomic phylogenetic trees from genbank files, it include automated closest reference detection, clustering along the references and generation of population level variant calling and variant annotation and generation of a tree of life like community level phylogenetic tree.",
    long_description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/ysevel/Clusttree",
    license='CCO 1.0 Universal',
    python_requires='=3.11',
    install_requires=[
                      'snippy>=4.6',
                      'bakta>=11.4',
                      "mash>=2.3.0",
                      "IQ-Tree>=2.4.0",
                      "getphylo>=1.0.1",
                      "Clonal_frame>=1.20",
                      "snp_dists>=1.2.0",
                      "ncbi_genome_download>=0.3.3"],
    entry_points={'console_scripts': ["Clusttree=src.phylogeny_workflow:main"]}
)
