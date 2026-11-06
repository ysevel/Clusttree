# Clusttree

Clusttree is a phylogenomic workflow for automatic reference selection, genome annotation, variant calling and phylogenetic reconstruction from bacterial genome assemblies.

The tool is designed for singlecell or metagenomic study. It will:

- determine the closest reference of each genome (SAG/MAG) using a very fast kmer method using mash

- cluster together all the genomes with the same reference, produce variant calling prediction on the reference and annotate each variant, it will finally produce a population level recombination aware phylogenetic tree using clonalframe and IQ-TREE

- produce a community level "tree of life" phylogeny using getphylo.

## Structure of the workflow:
<img width="1001" height="368" alt="Clustree_workflow" src="https://github.com/ysevel/Clusttree/Image/Workflow_structure.png](https://github.com/ysevel/Clusttree/blob/main/Image/Workflow_structure.png" />

## Parameters
### Required parameters

-i --input (required): Path to the folders containing the annotated genomes (format: genbank)

-o --output (required): path to the output folder.

-s --suffix (required): suffix for the output file and run directory.

-b --bakta_db (required): Path to bakta database for genome annotation

#### Mutually exclusive start parameters:

-r --reference_mesh_file : Reference sketch distance matrix file to use as a comparison with each SAG to determine the closest reference file.

-l --reference_list : File with the list of references to generate the sketch matrix from. The file should either be in the working directory or given as path.

-d --download-refseq : Output folder to download the whole refseq database for bacteria and archeae, will download all the reference genomes from refseq and automatically generate a sketch matrix for all those references.

### Optional parameters

-t --threads (default=1): number of threads allocated to the workflow.

-n --nb_closest_ref (default=1): Number of ref to be included. It will select the n closest reference from the query genome.

-ds --divide_step: Launch steps of the pipeline separately. By default launch all steps. Choices are "all", "reference_selection", "variant_calling" and "phylogenomy", see structure of the workflow for details about each step.



## Requirements:

Bakta >= 1.11.4: https://github.com/oschwengers/bakta

getphylo >= 1.0.1: https://github.com/drboothtj/getphylo

IQ-TREE >= 2.4.0: https://iqtree.github.io/

Mash >= 2.3.0: https://github.com/marbl/Mash

NCBI-genome-download >= 0.3.3: https://github.com/kblin/ncbi-genome-download

snippy >= 4.6.0: https://github.com/tseemann/snippy

ClonalFrameML >= 1.20: https://github.com/xavierdidelot/ClonalFrameML

snp-dists >= 1.20: https://github.com/tseemann/snp-dists 

