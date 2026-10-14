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

# Global variables
manager = Manager()
current_missing_id_list = manager.list() #keep as an option if needed to exclude SAgs from the analysis
missing_id_lock = Lock() #same as above, to remove if not included in the final worflow...

#functions
def handle_program_options():
    descr = """\
    Determination of the closest reference to each genome in a folder using mash distance from a reference genomes database. Also peform variant analysis and phylogenetic analysis on the genome based on their closest references.
    Program arguments:
    -----------------"""

    parser = argparse.ArgumentParser(description=twdd(descr),
                           formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help="path to the input folder containing the genbank files")
    parser.add_argument('-o', '--output', required=True,
                        help="Output path where directories will be built.")
    parser.add_argument('-t', '--threads', type=int ,default=1,
                        help="Number of threads")
    parser.add_argument('-s', '--suffix', required=True,
                        help="suffix added to each output file")
    parser.add_argument('-m', '--metadata_file', help="path to the reference database metadatas") #not used at the moment, obsolete? maybe turn into an output later on
    parser.add_argument('-n', '--nb_closest_ref', type=int, default=1, help="number of ref to be included. it will select the x closest reference from the query genome (default=1)") #keep to be able switch to mixed ref prediction
    parser.add_argument('-e', '--email', required=True, help='an email adress to be provided to Entrez as it is required to use E-utilities.') #maybe hardcode it in the code?
    group = parser.add_mutually_exclusive_group(required=True) #optional starting point for the pipeline, either build the db, or construct mash matrix or calculate the distance to an existing matrix.
    group.add_argument('-r', '--reference_mesh_file', help="reference mash mesh file to use as a comparison with each SAG to determine the closest reference file.")
    group.add_argument('-l', '--references_list', help= "file with the list of references to generate the mash mesh file from.")
    group.add_argument('-d', '--download_refseq', help='Output folder to download the whole refseq database for bacteria and archeae. Other databaese can be used, but have to be downloaded manually.')
    return parser.parse_args()

def generate_reference(path_ref_list, output_dir, suffix, nb_threads, logger):
    '''This function generate a mesh distance matrix from a list of reference genomes.'''
    logger.info("Running Mash distance sketch matrix")
    os.makedirs(output_dir, exist_ok=True)
    log_mash=f"{output_dir}/mash_dist_{suffix}.log"
    output_mesh_file=f"{output_dir}/NCBI_repr_{suffix}.msh" #output file get automatically named, it can be changed safely here.
    cmd_mash=f"mash sketch -p {nb_threads} -l {path_ref_list} -o {output_mesh_file}"
    subprocess.run(cmd_mash, shell=True)
    return output_mesh_file #return the mash sketch file to be used in launch mash.

def launch_mash(assembly_file, id_sample,ref_mash , output, suffix, threads, logger):
    """
    create the mash result directories and launch mash for each SAG on a reference mesh distance matrix.
    """
    os.makedirs(f"{output}/mash_distance", exist_ok=True)
    mash_dist_file=f"{output}/mash_distance/{id_sample}_dist_NCBI.txt" #format of the output individual distance matrix.
    cmd_mash=f"mash dist -p 4 {ref_mash} {assembly_file} > {output}/mash_distance/{id_sample}_dist_NCBI.txt"
    #logger.info(f"the mash command is: {cmd_mash}")
    log_mash=f"{output}/{id_sample}_{suffix}_filter_reads.log"
    results = subprocess.run(cmd_mash, capture_output=True, shell=True, check=True)
    #logger.info(results.stdout.strip()) #should look at the output in details in the log. is that the b'' in the log? to test before release...
    return mash_dist_file

def download_reference(ref_dir, suffix, path_ref_list, nb_threads, logger): 
    '''The function will automatically download the up to date refseq database for bacteria and archaea and run generate_reference on it to produce a new mesh distance matrix.'''  
    #download fasta and genbank files from NCBI refseq, designer to download the references from refseq, maybe add more custoim option later on
    os.makedirs(f"{ref_dir}", exist_ok=True)
    fasta_out_folder=f"{ref_dir}/fasta"
    logger.info(f"Downloading reference file from NCBI refseq, the fasta are stored in {fasta_out_folder}.")
    logger.info(f"Starting downloading refseq reference genome on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")
    cmd_ncbi_download_fasta=f"ncbi-genome-download bacteria,archaea -s refseq -F fasta -R reference -o {fasta_out_folder} --flat-output -r 2 -m {path_ref_list}/reference_db_{suffix}_metadata.tab"
    subprocess.run(cmd_ncbi_download_fasta, shell=True, check=True)
    logger.info(f"download of the referene sequences end. The metadata are available in {ref_dir}_{suffix}_refseq_metadata.tab")
    #unzip both fasta and genank files
    subprocess.run(f"gunzip -r {fasta_out_folder}", shell=True)
    #generate the reference file list into a text file
    with open(path_ref_list, 'w') as ref_path_f:
        fasta_files = [
        os.path.join(fasta_out_folder, f) for f in os.listdir(fasta_out_folder) 
        if f.endswith(('.fasta', '.fa', '.fna')) #use every fasta file in the fasta_out_folder to generate a reference list.
        ]
        for f in fasta_files:
            ref_path_f.write(f+"\n")
    reference_file=generate_reference(path_ref_list, output_dir, nb_threads, logger) #run generate_reference on the reference list.
    return reference_file #return the mesh distance matrix to be used as refence for the rest of the run.

def closest_ref_finder(id_sample, mash_dist_file, nb_closest_ref, reference_groups, output, logger):
    '''return the closest reference form a mash distance matrix, extract the accession for the reference annotation and metadata generation, note: it is recommanded to use references with the proper NCBI accession in their names.'''
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
        match = re.match(r'^(GCF|GCA)_\d+', file_name) #recover the accession ID from the file name.
        accession= None #accession do not exist by default.
    #Check if the accession exist, if so, save it and return it.
    if match:
        accession= match.group(0) #if a accession exist in the file name, return the accession.
        logger.info(f"The closest reference for {id_sample} is {accession}.") #save the accession in the log for each SAG, for verification steps.
    else:
        logger.error(f"Reference accession  not found for {id_sample}") #if no accession in the file name, return an error and write it in the log.
    return closest_ref_name, accession

def getphylo(genbank_dir,output_getphylo, seed, nb_threads, logger):
    '''Compute the community tree based on all the SAgs and all the predicted closest references using getphylo workflow, use the reference from the most abundant cluster as a seed to run the workflow.'''
    #run getphylo to produce a tree of life, get a selection of references to add to the tree.
    output_getphylo = os.path.abspath(output_getphylo)
    if output_getphylo.startswith(run_dir): #delete the folder if it already exist, getphylo does not work on previous run... check that output_getphylo start with run_dir for safety.
        shutil.rmtree(output_getphylo)
        logger.warning(f"The folder {output_getphylo} already exist, it has been deleted for getphylo to run again.")
    os.makedirs(f"{output_getphylo}", exist_ok=True)
    log_phylo=f"{output_getphylo}/getphylo_log_{suffix}.log"
    command_getphylo=f"getphylo -g '{genbank_dir}/*.gbff' -s {seed} -o {output_getphylo} -p 20 -c {nb_threads} -m iqtree -ia" #run get phylo on the unique gene present in 20% of the dataset.
    logger.info(f"Running getphylo with cmd: {command_getphylo}")
    with open(log_phylo, 'w') as log_p:
        subprocess.run(command_getphylo,stdout=log_p, stderr=subprocess.STDOUT, shell=True, check=True)
    return
    

def variantcalling(fasta_dir, output_dir, sag_list, nb_threads, closest_ref_gbk, suffix, logger): #run multi snippy on the dataset using the closest reference
    '''Function to peform variant calling and phylogeny on a cluster of strains using the same reference, the reference is used to align all assembly and perform variant calling, recombinasion free phylogeny is comptued using IQ-tree and ClonalFrameML '''
    log_snippy=f"{output_dir}/snippy_{suffix}.log"
    with open(log_snippy,'w') as log_v:
        logger.info(f"running snippy on each fasta file associated with {closest_ref_gbk}: {sag_list}")
        sag_id_list=[]
        reference_name=Path(closest_ref_gbk).stem
        #run snippy for each genome:
        for SAGs in sag_list:
            SAG_id=Path(SAGs).stem
            sag_id_list.append(SAG_id)
            command_snippy_variant=f"snippy --ctgs '{SAGs}' --ref {closest_ref_gbk} --cpus {nb_threads} --outdir {output_dir}/{SAG_id}"
            logger.info(f"running snippy on {SAG_id} with the command: {command_snippy_variant}")
            subprocess.run(command_snippy_variant,stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True)
            snpeff_ref_dir=f"{output_dir}/{SAG_id}/reference"
        #os.chdir(output_dir) #obsolete, can lead to various annoying bugs.
        sag_str=' '.join(f"{str(e)}" for e in sag_id_list) #check the output.
        command_snippy_core=f"snippy-core --ref {closest_ref_gbk} {sag_str} --prefix {suffix}_{reference_name}"
        logger.info(f"running snippy core for {reference_name}")
        logger.info(command_snippy_core)
        subprocess.run(command_snippy_core,cwd=output_dir, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check =True)
        logger.info(f"snippy core done, preparing subtree prediction")
        logger.info(f"compute mutation frequencies")
        #run IQ-tree to compute the initial tree
        logger.info(f"population genomic is done.")
        cmd_iqtree=f"iqtree -T AUTO -s {output_dir}/{suffix}_{reference_name}.full.aln -m MFP --prefix {suffix}_subtree_{reference_name} -safe" #use the maximal amount of threads.
        logger.info(f"running iqtree: {cmd_iqtree}")
        subprocess.run(cmd_iqtree, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True)
        #filter recombination event:
        logger.info("running ClonalFramel on the initial tree")
        cmd_clonalframe=f"ClonalFrameML {output_dir}/{suffix}_subtree_{reference_name}.treefile {output_dir}/{suffix}_{reference_name}.full.aln {output_dir}/{suffix}_{reference_name}_CFML -output_filtered true -ignore_incomplete_sites true" #produce filtered variant fasta
        logger.info(f"Clonalframe command: {cmd_clonalframe}")
        subprocess.run(cmd_clonalframe,cwd=output_dir , stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True) #ERROR: currently clonalframe won't run on genouest.
        logger.info("Running IQ tree on the recombination free alignement")
        #run IQ-tree on the filtered variants:
        cmd_iqtree_filtered=f"iqtree -T AUTO -s {output_dir}/{suffix}_{reference_name}_CFML.filtered.fasta -m MFP --prefix {suffix}_subtree_{reference_name} -safe"
        logger.info(f"final iq tree command: {cmd_iqtree_filtered}")
        subprocess.run(cmd_iqtree_filtered, stdout=log_v, stderr=subprocess.STDOUT, shell=True, check=True) #currently unavailable, must be put back on for the full intended workflow.
        logger.info(f"subtree succesfully computed for {closest_ref_gbk}")
        annotate_variants(output_dir, suffix, nb_threads, reference_name, snpeff_ref_dir, logger)
    return

def generate_chrom_map(config_file, output_map_file, logger):
    chromosomes = [] #creation of the chromosome list
    
    # find ref.chromosome and capture each chromosome in a list
    chrm_pattern = re.compile(r"^ref\.chromosome\s*:\s*(.*)$")

    try:
        with open(config_file, 'r') as f:
            for line in f:
                match = chrm_pattern.match(line.strip())
                if match:
                    # implementation of the chromosome list
                    raw_chroms = match.group(1).split(',')
                    chromosomes = [c.strip() for c in raw_chroms]
                    break
        
        if not chromosomes:
            logger.error("No ligne 'ref.chromosome' found in the config file")
            return

        # Generarion of the chromosome annotation file.
        with open(output_map_file, 'w') as out:
            for i, chrom_name in enumerate(chromosomes, start=1):
                out.write(f"{i} {chrom_name}\n")
        logger.info(f"correspondance table succesfully generated : {output_map_file}")

    except FileNotFoundError:
        logger.error(f"The contig file {config_file} cannot be found.")

def annotate_variants(path_subtree, suffix, nb_threads, reference_name, snpeff_ref_dir, logger):
    '''This function perform variant annotation based on the recombinsaion free fltered variants predicted by ClonalFrameML, the closest reference of each cluster of strain is used for the SNP annotation.'''

    #log config
    logger.info(f"variant annotation begin for {reference_name}.")
    log_annotation=f"{path_subtree}/VCF_annotation_{suffix}.log"

# process the annotations
    with open( log_annotation, 'w') as log_a:
        #generating the VCF file
        cmd_vcf=f"snp-sites -v -o {path_subtree}/{suffix}_{reference_name}_predicted_variants.vcf {path_subtree}/{suffix}_{reference_name}_CFML.filtered.fasta" #new temporary version
        #cmd_vcf=f"snp-sites -v -o {path_subtree}/{suffix}_{reference_name}_CFML_filtered_variants.vcf {path_subtree}/{suffix}_{reference_name}_CFML " #old version for Use with clonalframe replace the new with the old when the problem is fixed with clonal frame
        subprocess.run(cmd_vcf, stdout=log_a, stderr=subprocess.STDOUT, shell=True, check=True)
        snpeff_config_file= f"{snpeff_ref_dir}/snpeff.config"
        logger.info(f"snpeff config file: {snpeff_config_file}")
        chrom_map_file= f"{path_subtree}/renamme_chr.txt"
        generate_chrom_map(snpeff_config_file, chrom_map_file, logger)

        #generate the distance matrix
        cmd_ditance_matrix=f"snp-dists {path_subtree}/{suffix}_{reference_name}_CFML.filtered.fasta > {path_subtree}/{suffix}_{reference_name}_dist_matrix.tsv"
        subprocess.run(cmd_ditance_matrix, stdout=log_a, stderr=subprocess.STDOUT, shell=True, check=True)
        logger.info(f"Distance matrix generated: {path_subtree}/{suffix}_{reference_name}_dist_matrix.tsv")

        #annotate the SNP
        cmd_annotate=f"bcftools annotate --rename-chrs {chrom_map_file} {path_subtree}/{suffix}_{reference_name}_predicted_variants.vcf > {path_subtree}/{suffix}_{reference_name}_temp.vcf"
        subprocess.run(cmd_annotate, stdout=log_a, stderr=subprocess.STDOUT, shell=True, check=True)
        cmd_snpeff=f"snpEff -ud 0 -c {snpeff_config_file} -dataDir {snpeff_ref_dir} ref {path_subtree}/{suffix}_{reference_name}_temp.vcf > {path_subtree}/{suffix}_{reference_name}_annotated_variants.vcf"  
        subprocess.run(cmd_snpeff, stdout=log_a, stderr=subprocess.STDOUT, shell=True, check=True)
        cmd_add_allele_freq=f"bcftools +fill-tags {path_subtree}/{suffix}_{reference_name}_annotated_variants.vcf > {path_subtree}/{suffix}_{reference_name}_full_annotated_variants.vcf"
        subprocess.run(cmd_add_allele_freq, stdout=log_a, stderr=subprocess.STDOUT, shell=True, check=True)
        logger.info(f"variant annotation end for {reference_name}.")
    return
    
#commandes à intégrer


#execution du programme
def process_sample(path_input, suffix, ref_file, nb_ref, run_dir, nb_threads, email):
    ###initialization variables
    os.makedirs(f"{run_dir}", exist_ok=True)
    log_file = f"{run_dir}/closest_ref_predictor_{suffix}.log"
    print(f"All prints will now be saved in {log_file} file.")
    reference_list=[]
    SAGs_list=[]
    reference_groups = defaultdict(lambda: {'SAGs': [], 'accession': None, 'taxname': None, 'taxid': None, 'genbank_file': None})
    
    #Initialisation of the log
    if os.path.exists(log_file):    # generates a new file
        os.remove(log_file)
    #create and implement the log
    logger = logging.getLogger(suffix)
    logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    #Test the input path.
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
        ref_file=download_reference(path_ref, path_ref_list, suffix, nb_threads, logger)
    elif ref_list_path:
        path_ref_list=ref_list_path
        logger.info("The option generate ref mesh file is activated, running mash sketch on the references list")
        ref_file=generate_reference(path_ref_list, output_dir, nb_threads, logger)
    else:
        logger.info(f"Running mash on {ref_file}")

    #start processing samples
    fasta_folder=f"{run_dir}/fasta_files"
    genbank_output=f"{run_dir}/genbank_output"
    os.makedirs(f"{fasta_folder}", exist_ok=True)
    os.makedirs(f"{genbank_output}", exist_ok=True)
    #parse the genbnak files
    genbank_files = [
    os.path.join(path_input, f) for f in os.listdir(path_input) 
    if f.endswith(('.gb', '.gbk', '.gbff')) #hopefully, this will allow to find any fasta file in the input folder...
    ]
    summary_data= []
    fasta_files=[]
    logger.info(f"reference prediction finished. annotating closest refereneces")
    for gb_file in genbank_files:
        #convert genbank in fasta
        logger.info(f"Processing file: {gb_file}")
        shutil.copy(f"{gb_file}", f"{genbank_output}") #copy gb files to the genbank folder to avoid to put the references in the input folder.
        strain_name=Path(gb_file).stem
        fasta_file = f"{fasta_folder}/{strain_name}.fasta"
        #create fasta files for mash and snippy
        SeqIO.convert(gb_file, "genbank", fasta_file, "fasta")
        input_assembly= fasta_file
        fasta_files.append(fasta_file)
        id_sample=os.path.splitext(os.path.basename(fasta_file))[0]
        #reference prediction
        logger.info(f"assembly file found for {id_sample}, calculating mash distance\n")
        mash_dist_f=launch_mash(input_assembly, id_sample, ref_file, run_dir, suffix, nb_threads, logger)
        closest_reference,accession =closest_ref_finder(id_sample, mash_dist_f, nb_ref, reference_groups, run_dir, logger) #Keep only the closest reference and its accession
        reference_name=Path(closest_reference).stem
        #extract relevant metadata for each reference fro the accession determined by closest_ref_finder.
        taxname = subprocess.run(f"esearch -db assembly -query '{accession}' | esummary | xtract -pattern DocumentSummary -element SpeciesName", shell=True, capture_output=True, text=True).stdout.strip().replace(' ', '_')
        safe_taxname = re.sub(r'[^A-Za-z0-9._-]', '_', taxname) #generate a safe taxname to avoid strange characters that would mess up the workflow later on.
        taxid= subprocess.run(f"esearch -db assembly -query '{accession}' | esummary | xtract -pattern DocumentSummary -element Taxid", shell=True, capture_output=True, text=True).stdout.strip()
        summary_data.append({'id_sample': id_sample, 'reference': closest_reference, 'accession': accession, 'tax name': safe_taxname, 'taxid': taxid })
        #Add the closest reference to a list of reference, but will only do it once to avoid loosing some time.
        if closest_reference not in reference_list:
            reference_list.append(closest_reference) #Append the ref in the list only one time, in order to use them in the following steps of the pipeline
            logger.info(f"new reference added to the reference list: {closest_reference}, downloading the genbank file for {accession}.")
            annot_dir=f"{run_dir}/references_annotation"
            os.makedirs(f"{annot_dir}", exist_ok=True) #save the annotation for the reference in a separate folder.
            cmd_bakta=f"bakta -t {nb_threads} --output {annot_dir}/{accession}_annot -p {accession}_{safe_taxname} --force --compliant {closest_reference}" #Meant to use the default variable for bakta db... set it ahead of running the workflow...
            logger.info(f"running bakta: {cmd_bakta}")
            subprocess.run(cmd_bakta, shell=True, check=True) #generate the annotation for the reference, the file will be in the annotation folder named after the ref accession and taxname.
            shutil.copy(f"{annot_dir}/{accession}_annot/{accession}_{safe_taxname}.gbff", f"{genbank_output}") #copy each ref into the genbnak folder for getphylo
            gb_file=f"{annot_dir}/{accession}_annot/{accession}_{safe_taxname}.gbff"
            bakta_report=(f"{annot_dir}/{accession}_annot/{accession}_{safe_taxname}.txt")
            #implementation of the accession and taxname of a new reference to link the ref to its genbank file generated with bakta
            reference_groups[closest_reference]['accession']=accession
            reference_groups[closest_reference]['taxname']=safe_taxname
            reference_groups[closest_reference]['taxid']=taxid
            reference_groups[closest_reference]['genbank_file']=gb_file
    df_ref_summary=pd.DataFrame(summary_data)
    out_file_ref=f"{run_dir}/references_summary.csv" #save top ref into msh distance folder for tracability
    df_ref_summary.to_csv(out_file_ref, index=False)
    
    #creating the reference file output that give the number of SAGs for each references.
    references_groups=dict(reference_groups)
    ref_file=f"{run_dir}/references_list.csv"
    with open(ref_file , 'w') as ref_f:
        ref_f.write(f'reference,taxon name,nb_SAGs\n')
    print(f"reference annotation ended, processing variant calling with snippy on {len(reference_list)} subset.")
    logger.info(f"reference annotation ended, processing variant calling with snippy on {len(reference_list)} subset.")

    #this function will cluster the strains according to their closest reference, and run the variant calling function on each cluster.
    for reference in reference_list:
        ref_accession=reference_groups[reference]['accession']
        SAGs_list = reference_groups[reference]['SAGs']
        ref_gbk=reference_groups[reference]['genbank_file']
        ref_taxname=reference_groups[reference]['taxname']
        nb_SAGs=str(len(SAGs_list))
        reference_name=Path(reference).stem
        with open(ref_file , 'a') as ref_f:
            ref_f.write(f"{reference_name},{ref_taxname},{nb_SAGs}\n")
        output_snippy=f"{run_dir}/subtrees_{ref_accession}_{ref_taxname}" #output folder for the subtree.
        os.makedirs(f"{output_snippy}", exist_ok=True)
        #creation of the subtree for the current cluster.
        variantcalling(path_input, output_snippy, SAGs_list, nb_threads, ref_gbk, suffix, logger)
       
    #creating the tree of life using getphylo:
    #determine the reference that cluster the most SAGs.
    best_ref, best_ref_data = max(reference_groups.items(), key=lambda x: len(x[1]['SAGs']))
    best_reference_name=Path(best_ref).stem
    output_getphylo = f"{run_dir}/whole_collection_tree"
    os.makedirs(output_getphylo, exist_ok=True)
    seed=f"{genbank_output}/{best_ref_data['accession']}_{best_ref_data['taxname']}.gbff"
    logger.info(f"The most abundant reference is: {best_ref} ({len(best_ref_data['SAGs'])} SAGs)")
    logger.info(f"Running getphylo with {seed} as a seed")
    #run getphylo using the most abundant cluster as a seed.
    getphylo(genbank_output, output_getphylo, seed, nb_threads, logger)
    

#add a combination of the closest ref, cluster the strain according to their closest ref, return a list of closest ref and a list of strain associated to each ref.

if __name__ == "__main__": #note: à modifier pour supprimer checkM et passer le programme en souche par souche (input= fasta file)
    arguments = handle_program_options()
    path_input = arguments.input
    output = arguments.output           
    suffix = arguments.suffix
    ref_file=arguments.reference_mesh_file
    metadata=arguments.metadata_file
    nb_threads=arguments.threads
    nb_ref=arguments.nb_closest_ref
    email=arguments.email
    db_path=arguments.download_refseq
    ref_list_path=arguments.references_list

    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : closest reference selection start for {suffix}\n"
    print(message)    
    process_sample(path_input, suffix, ref_file, nb_ref, output, nb_threads, email)
    message = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : phylogeny end for {suffix}\n"
    print(message)

