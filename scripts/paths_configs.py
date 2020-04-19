#!/usr/bin/env python3
#DNAscan uses several tools and files. DNAscan paths to software binaries always endwith "/" while paths to files always end with the filename.
#Instructions about how to get the needed files are commented below near the coresponding path line.

import os.path

#Configs 

RG_ID = "4"

RG_LB = "lib1"

RG_PL = "illumina"

RG_PU = "unit1"

RG_SM = "20"

tmp_dir = ""

#this is the port used to reach the local deployment of the iobio serveces. E.g. if port = "8080" the iobio serveces will be available at http://localhost:8080 
port_num = "8080"

num_cpu = "1"


#Paths

path_to_db = ""

dnascan_dir = ""

path_gene_list = ""

path_iobio = "iobio/"

path_hisat = ""

path_samtools = ""

path_rtg = ""

path_multiqc = ""

path_fastqc = ""

path_sambamba = ""

path_freebayes = ""

path_annovar = ""

path_annovar_db = ""

#annovar_protocols and annovar_operations can be modified according to the wanted annovar databases. For instructions about
#how to download the annovar databases and to modify annovar_protocols and annovar_operations please visit 
#http://annovar.openbioinformatics.org/en/latest/ . 

annovar_protocols = "refGene,dbnsfp30a,clinvar_20170130,avsnp147"

annovar_operations = "g,f,f,f"

path_expansionHunter = ""

path_expansionHunter_jsons = ""

path_samblaster = ""

path_java = ""

path_vcftools = ""

path_scripts = "scripts/"

path_gatk = ""

path_bed = "data/test_data.bed"

#path_bed="/users/k1513213/brc_scratch/indels_project/gene_list_positions_sorted_no_overlap.bed"

#hg19 can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/
path_reference = ""

#hg19 index can be downloaded from ftp://ftp.ccb.jhu.edu/pub/data/hisat_indexes/hg19_hisat.tar.gz
#the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_hisat_index = ""

path_manta = ""

path_bedtools = ""

path_tabix = ""

#A viral index for the whole NCBI database of viral complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_index = ""

#A bacterial index for the whole NCBI database of bacterial complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_bacteria_index = ""

#The index can be created running "./$path_hisat/hisat-build $path_custom_reference index_base"
path_custom_microbes_index = ""

path_bwa = ""

path_bcftools = ""

#hg19 index can be downloaded from XXX
#the index can be created running "./$path_bwa/bwa index $path_reference"
path_bwa_index = ""

#custom tool options

hisat_custom_options = ""

bwa_custom_options = ""

freebayes_custom_options = ""

gatk_HC_custom_options = ""
