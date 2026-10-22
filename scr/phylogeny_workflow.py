##!/usr/bin/python
import os
from datetime import datetime
#import numpy as np #not used at the moment, keep as reference for potentiel improvements) 
import pandas as pd
import argparse
import subprocess
import shutil
import logging
import re
from collections import defaultdict
from pathlib import Path
from textwrap import dedent as twdd
from multiprocessing import cpu_count, Manager, Lock
from Bio import SeqIO
#from Bio import Entrez #I didn't used entrez cause it caused some bugs

# Global variables
manager = Manager()
current_missing_id_list = manager.list() #keep as an option if needed to exclude SAgs from the analysis
missing_id_lock = Lock() #same as above, to remove if not included in the final worflow...

#functions
def main():
    arguments = handle_program_options() #reads the argument parser
    
    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : closest reference selection start for {args.suffix}\n"
    print(message)
    
    path_input = arguments.input
    output = arguments.output           
    suffix = arguments.suffix
    ref_file=arguments.reference_mesh_file
    metadata=arguments.metadata_file
    nb_threads=arguments.threads
    nb_ref=arguments.nb_closest_refl
    db_path=arguments.download_refseq
    ref_list_path=arguments.references_list
   
    process_sample(path_input, suffix, ref_file, nb_ref, output, nb_threads, db_path=arguments.db_path, ref_list_path=arguments.ref_list_path)
    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : phylogeny end for {args.suffix}\n"
    print(message)

def handle_program_options():
    descr = """\
    Determination of the closest reference to each genome in a folder using mash distance from a reference genomes database. Also peform variant analysis and phylogenetic analysis on the genome based on their closest references.
    Program arguments:
    -----------------"""

    parser = argparse.ArgumentParser(description=twdd(descr),
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help="path to the input folder containing the fasta files")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built.")
    parser.add_argument('-t', '--threads', type=int ,default=1,
                        help="Number of threads")
    parser.add_argument('-s', '--suffix', required=True,
                        help="suffix added to each output file")
    parser.add_argument('-m', '--metadata_file', help="path to the reference database metadatas") #not used at the moment, obsolete? maybe turn into an output later on
    parser.add_argument('-n', '--nb_closest_ref', type=int, default=1, help="number of ref to be included. it will select the x closest reference from the query genome") #keep to be able switch to mixed ref prediction
    parser.add_argument('-g', '--genbank_dir', required=True, help='folder containing the genbank file associated with the assembly')
    parser.add_argument('-e', '--email', required=True, help='an email adress to be provided to Entrez as it is required to use E-utilities.') #maybe hardcode it in the code?
    group = parser.add_mutually_exclusive_group(required=True) #optional starting point for the pipeline, either build the db, or construct mash matrix or calculate the distance to an existing matrix.
    group.add_argument('-r', '--reference_mesh_file', help="reference mesh file to use as a comparison with each SAG to determine the closest reference file.")
    group.add_argument('-l', '--references_list', help= "file with the list of references to generate the mesh file from.")
    group.add_argument('-d', '--download_refseq', help='Output folder to download the whole refseq database for bacteria and archeae. The other section have to be downloaded manually.')
    return parser.parse_args()

def count_CDS(bakta_report, logger):
    out_file.write("genome\tCDSs\trRNA\thypothetical\n")
    for strain_id in list_strain:
        try:
            with open(bakta_report) as f:
                for line in f:
                    if("CDSs:" in line):
                        tar_cds=line.split(' ')[-1]
                        return tar_cds
        except:
            logger.warning(f"error processing bkata report for {ref_name}")
            return None

def generate_reference(path_ref_list, output_dir, nb_threads, logger):
    '''The function generate a mesh file from a list of reference genomes'''
    logger.info("Running Mash distance sketch matrix")
    os.makedirs(output_dir, exist_ok=True)
    log_mash=f"{output_dir}/mash_dist_{suffix}.log"
    output_mesh_file=f"{output_dir}/NCBI_repr_{suffix}.msh" #output file get automatically named, it can be changed safely here.
    cmd_mash=f"mash sketch -p {nb_threads} -l {path_ref_list} -o {output_mesh_file}"
    subprocess.run(cmd_mash, shell=True)
    return output_mesh_file #return the mash sketch file to be used in launch mash.

def launch_mash(assembly_file, id_sample,ref_mash , output, suffix, threads, logger):
    """
    create the result directories and launch mash for each SAG on a reference database.
    """
    os.makedirs(f"{output}/mash_distance", exist_ok=True)
    mash_dist_file=f"{output}/mash_distance/{id_sample}_dist_NCBI.txt"
    cmd_mash=f"mash dist -p 4 {ref_mash} {assembly_file} > {output}/mash_distance/{id_sample}_dist_NCBI.txt"
    #logger.info(f"the mash command is: {cmd_mash}")
    log_mash=f"{output}/{id_sample}_{suffix}_filter_reads.log"
    results = subprocess.run(cmd_mash, capture_output=True, shell=True)
    #logger.info(results.stdout.strip()) #should look at the output in details in the log. is that the b'' in the log? to test before release...
    return mash_dist_file

def download_reference(output_dir, suffix, path_ref_list, nb_threads, logger): 
    '''The function will automatically download the refseq database and update mash on it'''  
    #download fasta and genbank files from NCBI refseq, designer to download the references from refseq, maybe add more custoim option later on
    os.makedirs(f"{output_dir}", exist_ok=True)
    fasta_out_folder=f"{output_dir}/fasta"
    genbank_out_folder=f"{output_dir}/genbank"
    logger.info(f"Downloading reference file from NCBI refseq, the fasta are stored in {fasta_out_folder} and the genbank files in {genbank_out_folder}")
    logger.info(f"Starting downloading refseq reference genome on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")
    cmd_ncbi_download_fasta=f"ncbi-genome-download bacteria,archaea -s refseq -F fasta -R reference -o {fasta_out_folder} --flat-output -r 2 -m {path_ref_list}/reference_db_{suffix}_metadata.tab"
    #cmd_ncbi_download_genbank=f"ncbi-genome-download bacteria,archaea -s refseq -F genbank -R reference -o {genbank_out_folder} --flat-output -r 2 --species-taxids 562,543" #not used at this point
    subprocess.run(cmd_ncbi_download_fasta, shell=True)
    #subprocess.run(cmd_ncbi_download_genbank, shell=True) #same as above, not used at the moment
    logger.info(f"download of the referene sequences end. The metadata are available in {output_dir}_{suffix}_refseq_metadata.tab")
    #unzip both fasta and genank files
    subprocess.run(f"gunzip -r {fasta_out_folder}", shell=True)
    #subprocess.run(f"gunzip -r {genbank_out_folder}", shell=True)
    #generate the reference file list into a text file
    with open(path_ref_list, 'w') as ref_path_f:
        fasta_files = [
        os.path.join(fasta_out_folder, f) for f in os.listdir(fasta_out_folder) 
        if f.endswith(('.fasta', '.fa', '.fna'))
        ]
        for f in fasta_files:
            ref_path_f.write(f+"\n")
    reference_file=generate_reference(path_ref_list, output_dir, nb_threads, logger)
    return reference_file

def closest_ref_finder(id_sample, mash_dist_file, nb_closest_ref, reference_groups, output, logger):
    '''return the closest reference form a mash distance matrix'''
    with open(mash_dist_file, 'r') as mash_f:
        dist_mash= [line.strip().split("\t") for line in mash_f]
        df=pd.DataFrame( dist_mash, columns=["reference","query", "distance", "p-value", "shared_hashes"])
        df["distance"] = df["distance"].astype(float)
        top_ref=df.nsmallest(nb_closest_ref, "distance") #maybe not usefull any longer
        out_file=f"{output}/mash_distance/{id_sample}_top_ref.csv" #save top ref into msh distance folder for tracability
        top_ref.to_csv(out_file, index=False)
        logger.info(f"The {nb_closest_ref}th top ref are writen in the file: {out_file}") #do we keep it or not? can be used for fastANI...
        closest_ref=df.nsmallest(1, "distance").iloc[0]
        closest_ref_name=str(closest_ref["reference"])
        file_name=os.path.basename(closest_ref_name)
        sample_name=str(closest_ref["query"])
        reference_groups[closest_ref_name]['SAGs'].append(sample_name)
        match = re.match(r'^(GCF|GCA)_\d+', file_name) #look for accession id in order to recover the accession
        accession= None #more elegant way to deal with it, generate a problem with the reference retrieving later on.
    #Check if the accession exist, if so, save it and return it.
    if match:
        accession= match.group(0) 
        logger.info(f"The closest reference for {id_sample} is {accession}.") #save the accession in the log for each SAG, for verification steps.
    else:
        logger.error(f"Accession  not found for {id_sample}")
    return closest_ref_name, accession

def getphylo(genbank_dir,output_getphylo, seed, nb_threads, logger):
    #run getphylo to produce a tree of life, get a selection of references to add to the tree.
    if output_getphylo.startswith(run_dir): #delete the folder if it already exist, getphylo does not work on previous run... check that output_getphylo start with run_dir for safety. 
        shutil.rmtree(output_getphylo)
        logger.warning(f"The folder {output_getphylo} already exist, it has been deleted for getphylo to run again.")
    os.makedirs(f"{output_getphylo}", exist_ok=True)
    log_phylo=f"{output_getphylo}/getphylo_log_{suffix}.log"
    command_getphylo=f"getphylo -g '{genbank_dir}/*.gbff' -s {seed} -o {output_getphylo} -p 20 -c {nb_threads} -m iqtree -ia" #run get phylo on the unique gene present in 90% of the dataset.
    logger.info(f"Running getphylo with cmd: {command_getphylo}")
    with open(log_phylo, 'w') as log_p:
        subprocess.run(command_getphylo,stdout=log_p, stderr=subprocess.STDOUT, shell=True)
    return
    #this function will cluster the strains according to their closest reference, then output them to various list of strain, based on the reference.

def variantcalling(fasta_dir, output_dir, sag_list, nb_threads, closest_ref_gbk, suffix, logger): #run multi snippy on the dataset using the closest reference
    '''Function to peform variant calling and phylogeny on a cluster of strains using the same reference, the reference is used to align all assembly and perform variant calling, recombinasion free phylogeny is comptued using IQ-tree and ClonalFrameML '''
    log_snippy=f"{output_dir}/snippy_{suffix}.log"
    with open(log_snippy,'w') as log_v:
        logger.debug(f"running snippy on each fasta file associated with {closest_ref_gbk}: {sag_list}")
        sag_id_list=[]
        reference_name=Path(closest_ref_gbk).stem
        #run snippy for each genome:
        for SAGs in sag_list:
            SAG_id=Path(SAGs).stem
            sag_id_list.append(SAG_id)
            command_snippy_variant=f"snippy --ctgs '{SAGs}' --ref {closest_ref_gbk} --cpus {nb_threads} --outdir {output_dir}/{SAG_id}"
            logger.debug(f"running snippy on {SAG_id} with the command: {command_snippy_variant}")
            run_command(["snippy", "--ctgs", SAGs, "--ref", closest_ref_bgk, "--cpu", str(nb_threads), "--outdir", f"{output_dir}/{SAG_id}"])
            subprocess.run(command_snippy_variant,stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True)
            snpeff_ref_dir=f"{output_dir}/{SAG_id}/reference"
        #os.chdir(output_dir) #obsolete, can lead to various annoying bugs.
        sag_str=' '.join(f"{str(e)}" for e in sag_id_list) #check the output.
        command_snippy_core=f"snippy-core --ref {closest_ref_gbk} {sag_str} --prefix {suffix}_{reference_name}"
        logger.info(f"running snippy core for {reference_name}")
        logger.debug(command_snippy_core)
        subprocess.run(command_snippy_core,cwd=output_dir, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check =True)
        logger.info(f"snippy core done, preparing subtree prediction")
        cmd_iqtree=f"iqtree -T AUTO -s {output_dir}/{suffix}_{reference_name}.full.aln -m MFP --prefix {suffix}_subtree_{reference_name} -safe" #use the maximal amount of threads.
        logger.debug(f"running iqtree: {cmd_iqtree}")
        subprocess.run(cmd_iqtree, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True)
        #filter recombination event:
        logger.info("running ClonalFramel on the initial tree")
        cmd_clonalframe=f"ClonalFrameML {output_dir}/{suffix}_subtree_{reference_name}.treefile {output_dir}/{suffix}_{reference_name}.full.aln {output_dir}/{suffix}_{reference_name}_CFML -output_filtered true -ignore_incomplete_sites true" #produce filtered variant fasta
        logger.debug(f"Clonalframe command: {cmd_clonalframe}")
        subprocess.run(cmd_clonalframe,cwd=output_dir , stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True) #ERROR: currently clonalframe won't run on genouest.
        logger.info("Running IQ tree on the recombination free alignement")
        #run IQ-tree on the filtered variants:
        cmd_iqtree_filtered=f"iqtree -T AUTO -s {output_dir}/{suffix}_{reference_name}_CFML.filtered.fasta -m MFP --prefix {suffix}_subtree_{reference_name} -safe"
        logger.debug(f"final iq tree command: {cmd_iqtree_filtered}")
        subprocess.run(cmd_iqtree_filtered, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True) #currently unavailable, must be put back on for the full intended workflow.
        logger.info(f"subtree succesfully computed for {closest_ref_gbk}")
        annotate_variants(output_dir, suffix, nb_threads, reference_name, snpeff_ref_dir, logger)
    return

def process_sample(path_input, suffix, ref_file, nb_ref, run_dir, genbank_dir, nb_threads, email):
    ###initialize variable and directories
    os.makedirs(f"{run_dir}", exist_ok=True)
    log_file = f"{run_dir}/closest_ref_predictor_{suffix}.log"
    print(f"All prints will now be saved in {log_file} file.")
    reference_list=[]
    SAGs_list=[]
    reference_groups = defaultdict(lambda: {'SAGs': [], 'accession': None, 'taxname': None})
    
    ###initialization of the log
    if os.path.exists(log_file):    # generates a new file
        os.remove(log_file)
    #prepare log file
    logger = logging.getLogger(suffix)
    logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    if os.path.exists(path_input):
        logger.info(f"Running phylotools for {path_input} directory")
    else:
        logger.error("The sample directory do not exist, please check the path")
        return
    #check reference argument determine the starting point of the workflow.
    if db_path:
        path_ref= db_path
        logger.info(f"Downloading refseq database in {path_ref}")
        path_ref_list=f"{run_dir}/Path_ref_{suffix}_list.txt"
        ref_file=download_reference(path_ref, suffix, path_ref_list, nb_threads, logger)
    elif ref_list_path:
        path_ref_list=ref_list_path
        logger.info("The option generate ref mesh file is activated, running mash sketch on the references list")
        ref_file=generate_reference(path_ref_list, output_dir, nb_threads, logger)
    else:
        logger.info(f"Running mash on {ref_file}")

    #start processing samples
    fasta_files = [
    os.path.join(path_input, f) for f in os.listdir(path_input) 
    if f.endswith(('.fasta', '.fa', '.fna')) #hopefully, this will allow to find any fasta file in the input folder...
    ]
    for fasta_file in fasta_files:
        logger.info(f"Processing file: {fasta_file}")
        input_assembly= fasta_file
        id_sample=os.path.splitext(os.path.basename(fasta_file))[0]
        logger.info(f"assembly file found for {id_sample}, calculating mash distance\n")
        mash_dist_f=launch_mash(input_assembly, id_sample, ref_file, run_dir, suffix, nb_threads, logger)
        closest_reference,accession =closest_ref_finder(id_sample, mash_dist_f, nb_ref, reference_groups, run_dir, logger) #ne garder que la meilleurs reference...
        #add the closest reference to a list of reference, but will only do it once to avoid loosing some time.
        if closest_reference not in reference_list:
            reference_list.append(closest_reference) #append the ref in the list only one time, in order to use them in the following steps of the pipeline
            logger.info(f"new reference added to the reference list: {closest_reference}, downloading the genbank file for {accession}.")
            reference_name=Path(closest_reference).stem
            taxname = subprocess.run(f"esearch -db assembly -query '{accession}' | esummary | xtract -pattern DocumentSummary -element Organism", shell=True, capture_output=True, text=True).stdout.strip().replace(' ', '_')
            safe_taxname = re.sub(r'[^A-Za-z0-9._-]', '_', taxname)
            annot_dir=f"{run_dir}/references_annotation"
            os.makedirs(f"{annot_dir}", exist_ok=True) #save the annotation for the reference in a separate folder.
            cmd_bakta=f"bakta -t {nb_threads} --output {annot_dir}/{accession}_annot -p {accession}_{safe_taxname} --force --compliant {closest_reference}" #necessity to use the default variable for bakta db... set it ahead of running the workflow.
            logger.info(f"running bakta: {cmd_bakta}")
            subprocess.run(cmd_bakta, shell=True) #generate the annotation for the reference, the file will be in the annotation folder named after the ref accession and taxname.
            shutil.copy(f"{annot_dir}/{accession}_annot/{accession}_{safe_taxname}.gbff", f"{genbank_dir}") #copy each ref into the genbnak folder for getphylo
            bakta_report=(f"{annot_dir}/{accession}_annot/{accession}_{safe_taxname}.txt")
            #implementation of the accession and taxname of a new reference to link the ref to its genbank file generated with bakta
            reference_groups[closest_reference]['accession']=accession
            reference_groups[closest_reference]['taxname']=safe_taxname


    #creation of the subtree for each reference
    references_groups=dict(reference_groups)
    for reference in reference_list:
        SAGs_list=reference_groups.get(reference)
        reference_name=Path(reference).stem
        output_snippy=f"{run_dir}/subtrees_{reference_name}"
        os.makedirs(f"{output_snippy}", exist_ok=True)
        variantcalling(path_input, output_snippy, SAGs_list, nb_threads, reference, suffix, logger) 
        #ici la fonction qui va compter le nb de ref
       
    #phylgeny on genbank files for getphylo
    best_ref, best_ref_data = max(reference_groups.items(), key=lambda x: len(x[1]['SAGs']))
    best_reference_name=Path(best_ref).stem
    output_getphylo = f"{run_dir}/whole_collection_tree"
    os.makedirs(output_getphylo, exist_ok=True)
    seed=f"{genbank_dir}/{best_ref_data['accession']}_{best_ref_data['taxname']}.gbff"
    logger.info(f"The most abundant reference is: {best_ref} ({len(best_ref_data['SAGs'])} SAGs)")
    logger.info(f"Running getphylo with {seed} as a seed")
    getphylo(genbank_dir, output_getphylo, seed, nb_threads, logger)


if __name__ == "__main__": #initialisation of the workflow
    main()
