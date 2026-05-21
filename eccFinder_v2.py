
import os
import sys
import glob
import argparse
import multiprocessing
import subprocess

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from utilitaires.utilities import log,run_oae,get_eccFinder_version
from utilitaires.Aligner import Minimap2SAMAligner
from utilitaires.Aligner import Minimap2Aligner
from utilitaires.Peaker import genrich
from utilitaires.intersect import intersect_ecc,intersect_ecc_flexible

def read_genome_alignments(file_prefix,align_path,min_qlen,min_aln, overwrite_files):
    """ Filtering the raw read alignments based on query length and alignment length."""
    if os.path.isfile(align_path + file_prefix+".paf.bed"):
        log("INFO", "Retaining pre-existing file: " + align_path +file_prefix+".paf.bed")
    else:
        with open(align_path + file_prefix+".paf") as f:
            outfile=open(align_path + file_prefix+".paf.bed",'w')
            headers = ['refID','rstart', 'rend','queryID','qlen','direction']
            total=0
            n_primary=0
            for line in f:
                parts = line.strip().split("\t")
                total=total+1
                resultss = {
		        "queryID": parts[0],
		        "qlen":  int(parts[1]),
		        "qstart": int(parts[2]),
		        "qend": int(parts[3]),
		        "direction": parts[4],
		        "refID": parts[5],
		        "rlen": int(parts[6]),
		        "rstart": int(parts[7]),
		        "rend": int(parts[8]),
                "rdis": int(parts[8])-int(parts[7]),
		        "allmatch": int(parts[9]),
                "blockmatch": int(parts[10]),
                }
                if min_aln is not None and resultss['rdis'] < min_aln:
                    n_primary=n_primary+1
                    continue
                if min_qlen is not None and resultss['qlen'] < min_qlen:
                    continue
                out_row = (str(resultss[x]) for x in headers)
                outfile.write('\t'.join(out_row))
                outfile.write('\n')

def run_Porechop(file_prefix, query_file, peak_path, overwrite_files):
    """Run Porechop to trim adapters from reads."""
    query_fileporechop = os.path.join(peak_path, f"{file_prefix}.fastq")

    if os.path.isfile(query_fileporechop) and not overwrite_files:
        log("INFO", f"Retaining pre-existing file: {query_fileporechop}")
        return query_fileporechop

    log("INFO", f"Running Porechop on {query_file}")
    cmd = f"porechop -i {query_file} -o {query_fileporechop} --middle_threshold 75"
    subprocess.call(cmd, shell=True)

    return query_fileporechop


def run_TideHunter(file_prefix,query_fileporechop,peak_path, num_threads, max_divergence,min_period_size,num_copies,overwrite_files):
    #""" Spliting tandem repeats in one long read. """
    #query_fileporechop = os.path.join(peak_path, f"{file_prefix}.fastq")
    if os.path.isfile(peak_path +file_prefix+".unit.fa"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + peak_path +file_prefix+".unit.fa")
        else:
            log("INFO", "Overwriting pre-existing file: " +peak_path +file_prefix+".unit.fa")
            TH_params = " -c "+ str(num_copies)
            TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
            #TH_cmd = "./TideHunter/TideHunter-v1.5.5/bin/TideHunter"+ TH_params+ str(query_file)+ " -u > "+ peak_path +file_prefix+".unit.fa"

            TH_cmd = "TideHunter" + TH_params + str(query_fileporechop) + " -u > " + peak_path + file_prefix + ".unit.fa"

            subprocess.call(TH_cmd, shell=True)
    else:
        TH_params = " -c "+ str(num_copies)
        TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
        #TH_cmd = "./TideHunter/TideHunter-v1.5.5/bin/TideHunter"+ TH_params+ str(query_file)+ " -u > "+ peak_path +file_prefix+".unit.fa"
        TH_cmd = "TideHunter" + TH_params + str(query_fileporechop) + " -u > " + peak_path + file_prefix + ".unit.fa"
        subprocess.call(TH_cmd, shell=True)

def run_samtools(file_prefix, output_path, peak_path, num_threads, overwrite_files):
    """Sort, filter and index alignments properly for downstream peak calling."""

    final_bam = os.path.join(peak_path, f"{file_prefix}.unit.bam")
    input_sam = os.path.join(output_path, f"{file_prefix}.tmp.sam")
    modified_sam = os.path.join(output_path, f"{file_prefix}.illumina_like.sam")

    if os.path.isfile(final_bam) and not overwrite_files:
        log("INFO", f"Retaining pre-existing file: {final_bam}")
        return

    log("INFO", f"Generating BAM file: {final_bam}")

    # Step 0: Modify SAM to simulate Illumina-style read names
    log("INFO", f"Modifying SAM file: {input_sam} → {modified_sam}")
    with open(input_sam, "r") as infile, open(modified_sam, "w") as outfile:
        for line in infile:
            if not line.startswith("@"):
                fields = line.strip().split("\t")
                fields[0] = fields[0] + "/1"
                outfile.write("\t".join(fields) + "\n")
            else:
                outfile.write(line)

    # Step 1: Convert modified SAM to BAM
    tmp_bam = os.path.join(output_path, "tmp1.sam")
    pysam.view("-@", str(num_threads), "-bS", "-o", tmp_bam, modified_sam, catch_stdout=False)

    # Step 2: Sort BAM by query name
    sorted_bam = os.path.join(output_path, "tmp2.sam")
    pysam.sort("-@", str(num_threads), "-n", "-o", sorted_bam, tmp_bam, catch_stdout=False)

    # Step 3: Copy sorted BAM to final location
    subprocess.call(f"cp {sorted_bam} {final_bam}", shell=True)

    # Step 4: Index BAM
    #pysam.index(final_bam)


#######################################################################################################################################################

####################################################Function GeneRich ################################################################################

def run_Genrich(file_prefix,output_path,peak_path,num_threads, min_peak,max_dist,max_pvalue,overwrite_files):
    """ Detecting sites of genomic enrichment. """
    if os.path.isfile(output_path +file_prefix+".site.bed"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".site.bed")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".site.bed")
            run_samtools(file_prefix,output_path, num_threads, overwrite_files)
            GR_params = " -yv "
            GR_params += " -l " + str(min_peak)+" -g " + str(max_dist) +" -p " + str(max_pvalue)
            GR_cmd = "Genrich -t "+ peak_path +file_prefix+".unit.bam" + GR_params+ " -o "+ peak_path +file_prefix+".site"
            subprocess.call(GR_cmd, shell=True)
            cmd1 = "cut -f1-3 " + peak_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
            os.popen("{inS} ".format(inS=cmd1))
    else:
        GR_params = " -yv "
        GR_params += " -l " + str(min_peak)+" -g " + str(max_dist) +" -p " + str(max_pvalue)
        GR_cmd = "Genrich -t "+peak_path +file_prefix+".unit.bam" + GR_params+ " -o "+ peak_path +file_prefix+".site"
        subprocess.call(GR_cmd, shell=True)
        cmd1 = "cut -f1-3 " + peak_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
        os.popen("{inS} ".format(inS=cmd1))

##############################################################################################################################################################

##############################################FUNCTION TO FILTERED BED #####################################################################################################
def run_filterBED(file_prefix,output_path, align_path,min_read,min_bound,min_cov,overwrite_files):
    if os.path.isfile(output_path +file_prefix +".csv"):
        log("INFO", "Filtering locus by repeat units in one read")
    else:
        bedtools_params= " -wao -f "+ str(min_bound)
        #tmp1=output_path +file_prefix+".paf.bed.tmp1"
        tmp1=output_path +file_prefix+".paf.bed.txt"
        cmd ="bedtools intersect -a "+ output_path +file_prefix+".site.bed -b " + align_path +file_prefix + ".paf.bed" + bedtools_params+ " -nonamecheck > "+tmp1
        subprocess.call(cmd, shell=True)

        #cmd1 = "cat "+ output_path +file_prefix+".paf.bed.tmp1"
        cmd1 = "cat "+ output_path +file_prefix+".paf.bed.txt"
        cmd2 = "awk '$10>0'| sort -k7,7 -k9,9 | groupBy -g 1,2,3,7,9 -c 7 -o count |awk '$6>1'"
        cmd3 = "bedtools sort |groupBy -g 1,2,3 -c 4,6 -o count_distinct,sum -nonamecheck|awk '$4>2' > " +output_path +file_prefix+".paf.bed.tmp2"
        sub = "{inS} |{group}|{bed} ".format(inS=cmd1, group=cmd2, bed=cmd3)
        ps = subprocess.Popen(sub,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output = ps.communicate()[0]

        log("INFO", "Filtering locus by boundary coverage. ")
        csv=output_path +file_prefix +".tmp.csv"
        print(csv)
        cmd ="bedtools coverage -counts -a "+ output_path +file_prefix+".paf.bed.tmp2 -b " + align_path +file_prefix + ".paf.bed -f " + str(min_bound)+ " -nonamecheck > "+csv
        subprocess.call(cmd, shell=True)

        # Write the ecc_finder.csv.
        log("INFO", "Filtering locus by mininum number of tandely repeated long reads ")
        with open(output_path + file_prefix+".tmp.csv") as f:
            outfile=open(output_path + file_prefix+".csv",'w')
            headers = ['refID','rstart','rend','read_number','repeat_unit','coverage','ecc_len']
            total=0
            n_primary=0
            for line in f:
                parts = line.strip().split("\t")
                total=total+1
                resultss = {
	            "refID": parts[0],
	            "rstart":  int(parts[1]),
	            "rend": int(parts[2]),
	            "read_number": int(parts[3]),
	            "repeat_unit": parts[4],
                "coverage": int(parts[4]),
                "ecc_len":int(parts[2])-int(parts[1])
                }
                if min_read is not None and resultss['read_number'] < min_read:
                    n_primary=n_primary+1
                    continue
                if min_cov is not None and resultss['coverage'] < min_cov:
                    continue
                out_row = (str(resultss[x]) for x in headers)
                outfile.write('\t'.join(out_row))
                outfile.write('\n')

        for filename in glob.glob(os.path.join(output_path, "*tmp*")):
            try:
                #Trying to remove a current file
                os.remove(os.path.join(output_path, filename))
            except EnvironmentError:
                #You don't have permission to do it
                pass

#############################################################################################################################################################################

############################################Faire la Jointure des fichiers ###########################################################################
def joinFile(peak_path,file_prefix,output_path):
    df2 = pd.read_csv(peak_path +file_prefix + ".site", sep="\t", header=None)
    df1 = pd.read_csv(output_path + file_prefix+".csv",  sep="\t", header=None)
    print(df2[9])
    # Définir les colonnes de jointure par position
    colonnes_jointure = [0, 1, 2]  # positions des colonnes communes
    colonne_4_df2 = 7  # position de la 4e colonne à ajouter

    # Extraire les colonnes nécessaires du second fichier
    df2_reduit = df2[colonnes_jointure + [colonne_4_df2]]

    # Faire la jointure (INNER JOIN)
    fusion = pd.merge(df1, df2_reduit, on=colonnes_jointure, how="inner")
    fusion = fusion.drop(columns=[5])
    # 🔤 Ajouter tes propres noms de colonnes
    # Exemple : 3 colonnes clés + 2 colonnes supplémentaires
    noms_colonnes = ["chrom","chromStart","chromEnd", "Nb_reads", "Nb_hits","Peak_size","Pvalue"]  # adapte selon ton besoin
    fusion.columns = noms_colonnes

    # Sauvegarder le résultat
    fusion.to_csv(output_path +file_prefix +"_fusion.csv", sep="\t", index=False)
############################################################################################################################################################################
def runIntersect(file_prefix,output_path, gene_file,te_file,alu_file):
    intersect_ecc(
    ecc_file=output_path +file_prefix +"_fusion.csv",
    gene_annot=os.path.expanduser(gene_file),
    te_annot=os.path.expanduser(te_file),
    alu_annot=os.path.expanduser(alu_file),
    output_file=output_path +file_prefix +"final_output.tsv"
)

def runIntersect_flexible(file_prefix,output_path, gene_file):
    intersect_ecc_flexible(
    ecc_file=output_path +file_prefix +"_fusion.csv",
    annot_file=os.path.expanduser(gene_file),
    output_file=output_path +file_prefix +"final_output.tsv"
)


############################################################################################################################################################################
############################### FUNCTION GET RUN FASTA #####################################################################################################################
def run_getFasta(output_path, file_prefix ,ref_genome,overwrite_files):
    if os.path.isfile(output_path +file_prefix +".fasta"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".fasta")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".fasta")
            cmd ="seqtk subseq "+ ref_genome+ " " + output_path +file_prefix +".csv" + "> "+output_path +file_prefix +".fasta"
            subprocess.call(cmd, shell=True)
    else:
        cmd ="seqtk subseq "+ ref_genome+ " " + output_path +file_prefix +".csv" + "> "+output_path +file_prefix +".fasta"
        subprocess.call(cmd, shell=True)
###############################################################################################################################################################################

# function default thread number
def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)
# the parser is the function where we parse the argument to use especially function like tidehunter()
def main():
    #----- START TIDEHUNTER PARAMETER----- this part of argument it's necessary for run TideHunter function----------------------------------------------------------------
    description = "A tool to detect eccDNA loci using ONT sequencing"
    parser = argparse.ArgumentParser(description=description, usage="ecc_finder.py map-ont <reference.idx> or <reference.fai> <query.fq> -r <reference.fa> (option)")

    parser.add_argument("query", metavar="<query.fq>", nargs='?', default="", type=str,help="query fastq/fasta file (uncompressed or bgzipped)")

    parser.add_argument("idx", metavar="<reference.idx>", nargs='?', default="", type=str, help="index file of reference genome")
    parser.add_argument("-r", metavar="<query.fasta>", default="", type=str, help="reference genome fasta file (uncompressed or bgzipped)")
    #These 3 parameters below it's for gene, alu, te in gff file 
    parser.add_argument("-gene", metavar="STR", default="", type=str, help="annotated file for gene human or plant")
    parser.add_argument("-te", metavar="STR", default="", type=str, help="annotated file for te human or plant")
    parser.add_argument("-alu", metavar="STR", default="", type=str, help="annotated file for alu human or plant")
    parser.add_argument("-org", metavar="STR", default="human", type=str, help="organism  human or plant")
    val_options =parser.add_argument_group("validation options")
    val_options.add_argument("-n", metavar="INT", type=int, default=2, help="minimum copy number of tandem repeat in a long read [2]")
    val_options.add_argument("-e", metavar="FLT", type=float, default=0.25, help="maximum allowed divergence rate between two consecutive repeats [0.25]")
    val_options.add_argument("-s", metavar="INT", type=int, default=30, help="minimum period size of tandem repeat (>=2) [30]")    


     ################# PARAMETERS OF GENERICH ###############################################################
    peak_options = parser.add_argument_group("peak-calling options")
    peak_options.add_argument("-l", metavar="INT", type=int, default=200, help="minimum length of a peak [200]")
    peak_options.add_argument("-d", metavar="INT", type=int, default=100,
                              help="maximum distance between signif. sites [1000]")
    peak_options.add_argument("-p", metavar="FLT", type=float, default=0.05, help="maximum p-value [0.05]")
  

    #################################### PARAMETERS FILTERED BED #################################################
    val_options.add_argument("--min-read", metavar="INT", type=int, default=3, help="filter locus by unique mapped read number [3]")
    val_options.add_argument("--min-bound", metavar="FLT", type=float, default=0.8, help="filter locus at regions by boundary coverage (# aligned bases / boundary bases)[0.8]")
    val_options.add_argument("--min-cov", metavar="FLT", type=float, default=10, help="minimum coverage of detected eccDNA loci [10]")
    ################################################################################################################

    

    out_options = parser.add_argument_group("output options")
    out_options.add_argument("-o", metavar="PATH", type=str, default="eccFinder_output", help="output directory [./eccFinder_output]")
    out_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    out_options.add_argument("-x", type=str, default="ecc.ont", help="add prefix to output [ecc.ont]")
    out_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)
    map_options = parser.add_argument_group("map options")
    map_options.add_argument('-t', metavar="INT", type=int, default=get_default_thread(), help='number of CPU threads for mapping mode')
    mm2_default = "-x map-ont"
  
    #lire les arguments que l’utilisateur a fournis dans la ligne de commande et stocke-les dans args
    args = parser.parse_args()

    idx_file = args.r
    #default="ecc.ont", help="add prefix to output [ecc.ont]"
    file_prefix = args.x
    #help="query fastq/fasta file (uncompressed or bgzipped)
    query_file = os.path.abspath(args.query)
    # annotated file
    gene_file = os.path.abspath(args.gene)
    te_file = os.path.abspath(args.te)
    alu_file = os.path.abspath(args.alu)
    #default=2, help="minimum copy number of tandem repeat in a long read [2]"
    num_copies = args.n
    # default=0.25, help="maximum allowed divergence rate between two consecutive repeats [0.25]"
    max_divergence = args.e
    #default=30, help="minimum period size of tandem repeat (>=2) [30]"
    min_period_size = args.s
    #default=get_default_thread(), help='number of CPU threads for mapping mode'
    num_threads = args.t
    #default=False, help="overwrite intermediate files"
    overwrite_files = args.w
    ref_genome = args.r
    organism =args.org
    output_path = args.o
    # ----- END TIDEHUNTER PARAMETER---------------------------------------------------------------------------------------------------------------------------------
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"


    #default="eccFinder_output", help="output directory [./eccFinder_output]"
    peak_path = os.path.abspath(args.o) + "/peak_files/"
    if not os.path.isdir(output_path+"peak_files/"):
        os.makedirs(output_path+"peak_files/")
    peak_path = os.path.abspath(args.o) + "/peak_files/"
    #if not args.idx or not args.query or not args.r:
    if not args.query or not args.r:
        parser.print_help()
        sys.exit("\n** The reference fasta, idx and query files are required **")

    log("VERSION", "ecc_finder " + get_eccFinder_version())
    log("CMD", "python ecc_finder.py map-ont " + " ".join(sys.argv[1:]))
  
    ############################################FOLDERS FOR ALIGNMENTS#######################################################################
    align_path = os.path.abspath(args.o) + "/align_files/"
    if not os.path.isdir(output_path + "align_files/"):
        os.makedirs(output_path + "align_files/")
    align_path = os.path.abspath(args.o) + "/align_files/"
    ###################### PARAMETERS ALIGNMENTS#########################################"""
    map_options.add_argument("-q", metavar="INT", type=int, default=200, help="minimum query length [200]")
    map_options.add_argument("-a", metavar="INT", type=int, default=200, help="minimum alignment length [200]")
    args = parser.parse_args()
    min_qlen = args.q
    min_aln = args.a
    #############################################################################################################################################
    ################## Where you w  nt to launch alignments##################################################################
    # Align the query raw read to the reference.
    log("INFO", "Align the query raw read to the reference.")
    mm2_params = mm2_default
    mm2_params += " -t " + str(num_threads)
    mapont_aligner_path = "minimap2"
    log("INFO", "Mapping the query raw read to the reference genome")
    map_all = Minimap2Aligner(idx_file, [query_file], mapont_aligner_path, mm2_params, align_path + file_prefix,
                              in_overwrite=overwrite_files)
    print(map_all)
    map_all.run_aligner()

    # Filter raw read alignments based on query length and alignment length
    log("INFO", "Filtering read alignments based on query length and alignment length")
    read_genome_alignments(file_prefix, align_path, min_qlen, min_aln, overwrite_files)
     ##############################Launch TideHunter#############################################
    # Splitting into unit sequences of each tandem repeat for a long read
    log("INFO", "Detecting tandem repeat pattern from long reads")
    log("INFO", "./test_eccfinder.py")
    log("INFO", "Trimm adaptators")
    run_Porechop(file_prefix, query_file, peak_path, overwrite_files)
    log("INFO", "Trimmed adaptators")
    query_fileporechop = os.path.join(peak_path, f"{file_prefix}.fastq")
    run_TideHunter(file_prefix, query_fileporechop, peak_path, num_threads, max_divergence, min_period_size, num_copies,
                   overwrite_files)
    ######################################################################################
  
    # peak calling

    log("INFO", "Peak calling for reads with the tandem repeated pattern")
    minimap2_params = "-ax map-ont"
    mapont_aligner_path = "minimap2"
    mm2_params = minimap2_params + " -t " + str(num_threads)
    pre = output_path + file_prefix + ".tmp"
    map_all = Minimap2SAMAligner(idx_file, [peak_path + file_prefix + ".unit.fa"], mapont_aligner_path, mm2_params, pre,
                                 in_overwrite=overwrite_files)
    map_all.run_aligner()
    
    run_samtools(file_prefix,output_path,peak_path,num_threads, overwrite_files)
#run_Genrich(file_prefix, output_path, peak_path, num_threads, min_peak, max_dist, max_pvalue, overwrite_files)
    ################# PARAMETERS OF GENERICH ###############################################################
    #peak_options = parser.add_argument_group("peak-calling options")
    #peak_options.add_argument("-l", metavar="INT", type=int, default=200, help="minimum length of a peak [200]")
    #peak_options.add_argument("-d", metavar="INT", type=int, default=100,
    #                          help="maximum distance between signif. sites [1000]")
    #peak_options.add_argument("-p", metavar="FLT", type=float, default=0.05, help="maximum p-value [0.05]")
    # appell args pour ajouter des options
    #args = parser.parse_args()
    #############################
    min_peak = args.l
    max_dist = args.d
    max_pvalue = args.p
    ######### call GeneRich##################################################################################
    run_Genrich(file_prefix, output_path, peak_path, num_threads, min_peak, max_dist, max_pvalue, overwrite_files)

    #################################### PARAMETERS FILTERED BED #################################################
    #val_options.add_argument("--min-read", metavar="INT", type=int, default=3, help="filter locus by unique mapped read number [3]")
    #val_options.add_argument("--min-bound", metavar="FLT", type=float, default=0.8, help="filter locus at regions by boundary coverage (# aligned bases / boundary bases)[0.8]")
    #val_options.add_argument("--min-cov", metavar="FLT", type=float, default=10, help="minimum coverage of detected eccDNA loci [10]")
    #########appel args pour ajouter des options pour trié les fichiers Bed######################################
    #args = parser.parse_args()
    #############################################"
    min_read = args.min_read
    min_bound = args.min_bound
    min_cov = args.min_cov

    if min_read < 0:
        if min_read != -1:
            raise ValueError("--min-read must be >=3")
    if min_bound < 0:
        if min_bound != -1:
            raise ValueError("--min-cov must be >=0")
    if ref_genome:
        ref_genome = os.path.abspath(ref_genome)

    # Peak calling
    log("INFO", "Producing bed file of eccDNA locus")
    run_filterBED(file_prefix, output_path, align_path, min_read, min_bound, min_cov, overwrite_files)
    log("INFO", "End of prod  bed file of eccDNA locus")
#run_filterBED(file_prefix, output_path, align_path, min_read, min_bound, min_cov, overwrite_files)
    log("INFO", "JOIN FILE file of eccDNA locus")
    joinFile(peak_path, file_prefix, output_path)
    log("INFO", "INTERSECT FILE file of eccDNA")
    if organism=="human":
    	runIntersect(file_prefix,output_path,gene_file,te_file,alu_file)
    if organism=="plant":
        runIntersect_flexible(file_prefix,output_path,gene_file)
    log("INFO", "Producing fasta file of eccDNA locus")
    run_getFasta(output_path, file_prefix, ref_genome, overwrite_files)



if __name__ == '__main__':
    main()
