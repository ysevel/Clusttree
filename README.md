# Singelcell_phylogenic_workflow

Clusttree is a phylogenomic workflow for automatic reference selection, genome annotation, variant calling and phylogenetic reconstruction from bacterial genome assemblies. it's currently a very beta version that will undergo various updates soon.

The tools is designed study a diverse community from singlecell or metagenomic study. it will:

- determine the closest reference of each genome (SAG/MAG) using a very fast kmer method using mash)

- cluster together all the genome with the same reference, produce variant calling prediction on the reference and annotate each variant, it will finally produce a population level recombination aware phylogenetic tree. using clonalframe and IQ tree

- it will produce a community level "tree of life" phylogeny using getphylo.

# Structure of the workflow:
<img width="1001" height="368" alt="Clustree_workflow" src="https://github.com/user-attachments/assets/6eded873-f917-4978-b014-1ac9693fa6d0" />
