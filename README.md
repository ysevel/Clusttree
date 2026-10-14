# Singelcell_phylogenic_workflow
this is a beta version of a workflow for phylogeneic reconstruction.

This tool aim at reconstruct the population structure and the phylogeny of a bacterial community.

The tools is designed study a diverse community from singlecell of metagenomic study. it will:

- determine the closest reference of each genome (SAG/MAG) using a very fast kmer method using mash)

- cluster together all the genome with the same reference, produce variant calling prediction on the reference and annotate each variant, it will finally produce a population level recombination aware phylogenetic tree. using clonalframe and IQ tree

- it will produce a community level "tree of life" phylogeny using getphylo.

# Structure of the workflow:
<img width="1001" height="368" alt="Clustree_workflow" src="https://github.com/user-attachments/assets/6eded873-f917-4978-b014-1ac9693fa6d0" />
