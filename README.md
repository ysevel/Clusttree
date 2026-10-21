# Singelcell_phylogenic_workflow

Clusttree is a phylogenomic workflow for automatic reference selection, genome annotation, variant calling and phylogenetic reconstruction from bacterial genome assemblies. it's currently a very beta version that will undergo various updates soon.

The tools is designed study a diverse community from singlecell or metagenomic study. it will:

- determine the closest reference of each genome (SAG/MAG) using a very fast kmer method using mash)

- cluster together all the genome with the same reference, produce variant calling prediction on the reference and annotate each variant, it will finally produce a population level recombination aware phylogenetic tree. using clonalframe and IQ tree

- it will produce a community level "tree of life" phylogeny using getphylo.

# Structure of the workflow:
<img width="1001" height="368" alt="Clustree_workflow" src="https://github.com/user-attachments/assets/6eded873-f917-4978-b014-1ac9693fa6d0" />

# parameters

-i --input (required): Path to the folders containing the annotated genomes (format: genbank)
-o --output (required): path to the output folder.
-s --suffix (required): suffix for the output file and run directory.
-e email (required): an email address to be provided to Entrez as it is required to use E-utilities.
-t –threads (default=1): number of threads allocated to the workflow.
Mutually exclusive start parameters:
-r –reference_mesh_file : Reference sketch distance matrix file to use as a comparison with each SAG to determine the closest reference file.
-l –reference_list : File with the list of references to generate the sketch matrix from. The file should either be in the working directory or given as path.
-d –download-refseq : Output folder to download the whole refseq database for bacteria and archeae, will download all the reference genomes from refseq and automatically generate a sketch matrix for all those references.

# requierments:

Bakta: 1.11.4 https://github.com/oschwengers/bakta
Diamond: 2.1.10 https://github.com/bbuchfink/diamond
FastTree: 2.2.0 https://github.com/morgannprice/fasttree
getphylo: 1.0.1 https://github.com/drboothtj/getphylo
IQ-Tree: 2.4.0 https://iqtree.github.io/
Mash: 2.3.0 https://github.com/marbl/Mash
Muscle: 3.8.1551 https://github.com/rcedgar/muscle
NCBI-genome-download: 0.3.3 https://github.com/kblin/ncbi-genome-download
snippy: 4.6.0 https://github.com/tseemann/snippy
ClonalFrameML: 


