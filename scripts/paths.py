#!/home/aga/ENV/bin/python
#DNAscan uses several tools and files. DNAscan paths to software binaries always endwith "/" while paths to files always end with the filename.
#Instructions about how to get the needed files are commented below near the coresponding path line.

import os.path

path_iobio = "iobio/"

path_hisat = "/users/k1513213/brc_scratch/annotation_unility_project/hisat2-2.1.0/"

path_samtools = ""

path_rtg = ""

path_multiqc = ""

path_fastqc = ""

path_freebayes = "/opt/apps/bioinformatics/freeBayes/v1.0.2-33-gdbb6160-dirty/bin/"

path_annovar = "~/brc_scratch/local/annovar/"

path_annovar_db = "~/brc_scratch/local/annovar/humandb/"

path_expansionHunter = ""

path_expansionHunter_jsons = "/users/k1513213/brc_scratch/annotation_unility_project/json_2/"


#this is the port used to reach the local deployment of the iobio serveces. E.g. if port = "8080" the iobio serveces will be available at http://localhost:8080 
port_num = "8080"

path_java = ""

path_vcftools = ""

path_gatk = "/opt/apps/bioinformatics/gatk/3.5/"

num_cpu = "8"

path_bed = "/users/k1513213/brc_scratch/annotation_unility_project/WGS_Miseq_test/miseq/miseq_nochr.bed"

#path_bed="/users/k1513213/brc_scratch/indels_project/gene_list_positions_sorted_no_overlap.bed"


#hg19 can be downloaded from 
path_reference = "/mnt/lustre/datasets/projectmine/reference/Homo_sapiens/reference_data/Homo_sapiens/GRCh37-lite/GRCh37-lite.fa"

#hg19 index can be downloaded from ftp://ftp.ccb.jhu.edu/pub/data/hisat_indexes/hg19_hisat.tar.gz
#the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_hisat_index = "/mnt/lustre/datasets/projectmine/reference/Homo_sapiens/reference_data/Homo_sapiens/GRCh37-lite/GRCh37-lite"

path_manta = "/users/k1513213/brc_scratch/manta/manta-1.1.1.centos5_x86_64/bin/"

path_bedtools = ""

path_tabix = ""

#A viral index for the whole NCBI database of viral complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_index = "/users/k1513213/brc_scratch/annotation_unility_project/hisat2-2.1.0/index/viral.genome"

#A bacterial index for the whole NCBI database of bacterial complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_bacteria_index = ""

path_bwa = ""

path_bcftools = ""

#hg19 index can be downloaded from XXX
#the index can be created running "./$path_bwa/bwa index $path_reference"
path_bwa_index = "/mnt/lustre/datasets/projectmine/reference/Homo_sapiens/reference_data/Homo_sapiens/GRCh37-lite/GRCh37-lite.fa"
