#!/usr/bin/env python3

################################################################
# Program: DNAscan
# Version 0.1
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

################################################################
# Script structure:
# 1. Import python modules
# 2. Define paths_configs viriables
# 3. Define options from command line
# 4. Parse options from command line
# 5. Create working dir tree
# 6. Bed splitting
# 7. Remove duplicates
# 8. Alignment
#   8.1 Aligns paired end reads
#       8.1.1 Fast mode alignment
#       8.1.2 Normal and intesive mode alignment
#   8.2 Aligns single end reads
#       8.2.1 Fast mode alignment
#       8.2.2 Normal and intesive mode alignment
# 9. If input file is a sam file, it converts it into bam
# 10. Variant (snv and indel) calling
#   10.1 GATK hc indel calling (only performed in intensive mode)
#       10.1.1 identification of potential indel sites
#       10.1.2 GATK hc indel calling on selected positions (from previous step 10.1.1)
#   10.2 Freebayes snv and indel calling
# 11. Perform variant hard filtering
# 12. Perform known expansions search with ExpansionHunter
# 13. Structural Variant calling
# 14. Annotation with Annovar
# 15. Microbes screening
#   15.1 Exctract non human reads
#   15.2 Identifies present non-human microbes
#       15.2.1 Identifies present viruses
#       15.2.2 Generates virus report
#       15.2.3 Identifies present bacteria
#       15.2.4 Generates bacteria report
#       15.2.5 Identifies present user-selected microbes
#       15.2.6 Generates user-selected microbes report
# 16. Alignment report generation ( samtools flagstat and stats )
# 17. Sequencing data report generation ( fastqc )
# 18. Snv and indel calling report generation ( bcftools stats )
# 19. Html report generation ( Multiqc )
# 20. Annotated variants report generation
# 21. Starting iobio services
#################################################################

# 1. Import necessary python modules

import argparse
import sys
import string
import os
import os.path
import re
import paths_configs
import subprocess

from argparse import RawTextHelpFormatter

# 2. Define paths_configs variables from paths_configs.py

path_iobio = paths_configs.path_iobio

path_gene_list = paths_configs.path_gene_list

port_num = paths_configs.port_num

path_vcftools = paths_configs.path_vcftools

path_gatk = paths_configs.path_gatk

gatk_HC_custom_options = paths_configs.gatk_HC_custom_options

path_multiqc = paths_configs.path_multiqc

path_java = paths_configs.path_java

path_fastqc = paths_configs.path_fastqc

path_hisat = paths_configs.path_hisat

path_hisat_index = paths_configs.path_hisat_index

hisat_custom_options = paths_configs.hisat_custom_options

path_samtools = paths_configs.path_samtools

path_freebayes = paths_configs.path_freebayes

freebayes_custom_options = paths_configs.freebayes_custom_options

path_annovar = paths_configs.path_annovar

path_annovar_db = paths_configs.path_annovar_db

annovar_protocols = paths_configs.annovar_protocols

annovar_operations = paths_configs.annovar_operations

num_cpu = paths_configs.num_cpu

path_bed = paths_configs.path_bed

path_to_db = paths_configs.path_to_db

path_expansionHunter = paths_configs.path_expansionHunter

path_reference = paths_configs.path_reference

path_expansionHunter_jsons = paths_configs.path_expansionHunter_jsons

path_manta = paths_configs.path_manta

path_bedtools = paths_configs.path_bedtools

path_tabix = paths_configs.path_tabix

path_virus_index = paths_configs.path_virus_index

path_bacteria_index = paths_configs.path_virus_index

path_custom_microbes_index = paths_configs.path_custom_microbes_index

path_bwa = paths_configs.path_bwa

path_bwa_index = paths_configs.path_bwa_index

bwa_custom_options = paths_configs.bwa_custom_options

path_rtg = paths_configs.path_rtg

dnascan_dir = paths_configs.dnascan_dir

path_sambamba = paths_configs.path_sambamba

path_bcftools = paths_configs.path_bcftools

RG_ID = paths_configs.RG_ID

RG_LB = paths_configs.RG_LB

RG_PL = paths_configs.RG_PL

RG_PU = paths_configs.RG_PU

RG_SM = paths_configs.RG_SM

path_scripts = paths_configs.path_scripts

path_samblaster = paths_configs.path_samblaster

# 3. Define options variables from command line

parser = argparse.ArgumentParser(
    prog='python DNAscan.py',
    usage=
    '%(prog)s [options] -format "string" -reference "string" -in "string" \n please use the --help flag to print further information\n ',
    description=
    '############DNAscan help Message############ \n\nDNAscan uses the file paths_configs.py to locate the needed tools and files. Please make sure your paths_configs.py file is properly filled \n\nUsage example: \n\npython DNAscan.py -format fastq -out /home/user/test/ -in sample_sorted_1.fq.gz -reference hg19 -alignment -variantcalling -results_report\n\nPlease check the following list of optional and required options\n\n################################################',
    formatter_class=RawTextHelpFormatter)

parser.add_argument(
    '-RG',
    action="store",
    dest="RG",
    default=False,
    help=
    'if this flag is set the alignment stage will use the read group provided in paths_configs.py (Default = "False")'
)

parser.add_argument(
    '-alsgenescanner',
    action="store_true",
    dest="alsgenescanner",
    default=False,
    help=
    'if this flag is set DNAscan will perform the ALSgeneScanner analyses and reports (Default = "False")'
)

parser.add_argument(
    '-mode',
    action="store",
    dest="mode",
    default="fast",
    help='options are fast, normal, intensive [string] (default = "fast")')

parser.add_argument(
    '-filter_string',
    action="store",
    dest="filter_string",
    default="",
    help='bcftools filter string, eg "GQ>20 & DP>10" (Default = "")')

parser.add_argument(
    '-paired',
    action="store",
    dest="paired",
    default="1",
    help=
    'options are 1 for paired end reads and 0 for single end reads (Default = "1")'
)

parser.add_argument(
    '-vcf', action="store", dest="vcf", help='complementary vcf file')

parser.add_argument(
    '-in2',
    action="store",
    dest="input_file2",
    help='input file 2, for paired end reads only (fastq file)')

parser.add_argument(
    '-iobio',
    action="store_true",
    dest="iobio",
    help=
    'if this flag is set the iobio services will be started at the end of the analysis (Default = "False")',
    default=False)

parser.add_argument(
    '-alignment',
    action="store_true",
    dest="alignment",
    help=
    'if this flag is set the alignment stage will be performed (Default = "False")',
    default=False)

parser.add_argument(
    '-expansion',
    action="store_true",
    dest="expansion",
    help=
    'if this flag is set DNAscan will look for the expansions described in the jason folder in paths_configs.py  (Default = "False") ',
    default=False)

parser.add_argument(
    '-SV',
    action="store_true",
    dest="SV",
    help=
    'if this flag is set the structural variant calling stage will be performed (Default = "False") ',
    default=False)

parser.add_argument(
    '-BED',
    action="store_true",
    dest="BED",
    help=
    'restrict the analysis to the regions in the bed file (Default = "False") ',
    default=False)

parser.add_argument(
    '-virus',
    action="store_true",
    dest="virus",
    help=
    'if this flag is set DNAscan will perform viral scanning (Default = "False")  ',
    default=False)

parser.add_argument(
    '-bacteria',
    action="store_true",
    dest="bacteria",
    help=
    'if this flag is set DNAscan will perform bacteria scanning (Default = "False") ',
    default=False)

parser.add_argument(
    '-custom_microbes',
    action="store_true",
    dest="custom_microbes",
    help=
    'if this flag is set DNAscan will perform a customized microbe scanning according to the provided microbe data base in paths_configs.py (Default = "False")  ',
    default=False)

parser.add_argument(
    '-variantcalling',
    action="store_true",
    dest="variantcalling",
    help=
    'if this flag is set DNAscan will perform snv and indel calling (Default = "False")  ',
    default=False)

parser.add_argument(
    '-annotation',
    action="store_true",
    dest="annotation",
    help=
    'if this flag is set DNAscan will annotate the found variants (Default = "False")  ',
    default=False)

parser.add_argument(
    '-results_report',
    action="store_true",
    dest="results_report",
    help=
    'if this flag is set DNAscan generate a results report (Default = "False")  ',
    default=False)

parser.add_argument(
    '-alignment_report',
    action="store_true",
    dest="alignment_report",
    help=
    'if this flag is set DNAscan generate an alignment report (Default = "False") ',
    default=False)

parser.add_argument(
    '-sequencing_report',
    action="store_true",
    dest="sequencing_report",
    help=
    'if this flag is set DNAscan generate a report describing the input sequencing data (Default = "False") ',
    default=False)

parser.add_argument(
    '-calls_report',
    action="store_true",
    dest="calls_report",
    help=
    'if this flag is set DNAscan generate a report describing the found snvs and small indels (Default = "False")',
    default=False)

parser.add_argument(
    '-sample_name',
    action="store",
    dest="sample_name",
    default="sample",
    help='specify sample name [string] (default = "sample")')

parser.add_argument(
    '-rm_dup',
    action="store",
    dest="rm_dup",
    help=
    'do you want to mark duplicates while aligning the reads to the reference genome? Options are True or False [string] (Default = "True")',
    default="True")

parser.add_argument(
    '-debug',
    action="store_true",
    dest="debug",
    help=
    'if this flag is set DNAscan will not delete intermediete and temporary files (Default = "False")',
    default=False)

parser.add_argument(
    '-exome',
    action="store_true",
    dest="exome",
    help=
    'if this flag is set DNAscan will only look for variants in the whole exome (Default = "False")',
    default=False)

parser.add_argument(
    '-format',
    action="store",
    dest="format",
    default="fastq",
    help='options are bam, sam, fastq, vcf [string] ')

parser.add_argument(
    '-reference',
    action="store",
    dest="reference",
    default="hg19",
    help=
    'options are hg19, hg38, grch37 and grch38 the path to the reference fasta file must be specified in paths_configs.py [string]'
)

parser.add_argument(
    '-ref_file',
    action="store",
    dest="ref_file",
    default= "",
    help=
    'path to the reference file in fasta format [string]'
)

parser.add_argument(
    '-dnascan_dir',
    action="store",
    dest="dnascan_main_dir",
    default= "",
    help=
    'path to the DNAscan main dir [string]'
)





requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument(
    '-out',
    required=True,
    action="store",
    dest="out",
    help=
    'path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/'
)


requiredNamed.add_argument(
    '-in',
    required=True,
    action="store",
    dest="input_file",
    help='input file [string]')

# 4. Parse command line options

args = parser.parse_args()

mode = args.mode

alsgenescanner = args.alsgenescanner

exome = args.exome

format = args.format

paired = args.paired

reference = args.reference

input_file = args.input_file

input_file2 = args.input_file2

expansion = args.expansion

out = args.out + '/'

vcf = args.vcf

SV = args.SV

BED = args.BED

bacteria = args.bacteria

virus = args.virus

alignment = args.alignment

iobio = args.iobio

variantcalling = args.variantcalling

annotation = args.annotation

results_report = args.results_report

alignment_report = args.alignment_report

sequencing_report = args.sequencing_report

calls_report = args.calls_report

sample_name = args.sample_name

rm_dup = args.rm_dup

debug = args.debug

filter_string = args.filter_string

custom_microbes = args.custom_microbes

RG = args.RG

ref_file = args.ref_file 

dnascan_main_dir = args.dnascan_dir

if ref_file:
    
    path_reference = ref_file
    
if dnascan_main_dir:
    
    dnascan_dir = dnascan_main_dir

variant_results_file = ""

if format == "vcf":

    variant_results_file = input_file

    if not vcf:

        vcf = variant_results_file
        

print("\n################################################")
        
        

# 5. Create working dir tree

os.system(
    "mkdir %s ; mkdir  %slogs ; mkdir  %sreports ; mkdir  %sresults ; mkdir %stmp"
    % (out, out, out, out, out))


5.1

print(
                "############DNAscan Options############ \n\n DNAscan is running an anlysis with the following specifics:\n"
            )

options_log = open('%s/logs/options.log' %(out), 'w')

for arg in vars(args):

        print(arg,':    ', getattr(args,arg))
        
        options_log.write(arg + ':    ' + str(getattr(args,arg)) + "\n")
   
options_log.close()

print('\nOptions saved onto %s/logs/options.log \n' %(out))
      
      
# 6. Bed splitting: splitting the analysis region into subsets of equal length to distribute the work across the available threads.
# To do this DNAscan uses a bed file.
# If bed file is not provided, it generates one starting from the
# reference genome index. The pipeline uses a bed file to split the
# analyses in several subprocesses

# Y. adapt DB to reference

if reference == "grch37" or  reference == "grch38" :
    
    if reference == "grch37":
        
        ref_hg = "hg19"
        
    else:
        
        ref_hg = "hg38"
    
    os.system("zcat %s/exome_%s.bed.gz | sed 's/chr//g' | bgzip -c > %s/exome_%s.bed.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("zcat %s/%s_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > %s/%s_gene_db.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("cp %s/%s_gene_names.txt.gz  %s/%s_gene_names.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    
if alsgenescanner:

    alignment = True

    variantcalling = True

    annotation = True

    #SV = True

    #expansion = True

    BED = True

    path_expansionHunter_jsons = '%s/repeats/%s/' % (dnascan_dir, reference)

    path_bed = '%s/als_gene_scanner.bed' % (path_to_db)

    path_gene_list = ""

    annovar_operations = "g,f,f"

    annovar_protocols = "refGene,dbnsfp30a,clinvar_20170905,"
    
# Y. adapt DB to reference

if reference == "grch37" or  reference == "grch38" :
    
    
    
    if reference == "grch37":
        
        ref_hg = "hg19"
        
    else:
        
        ref_hg = "hg38"
    
    os.system("zcat %s/exome_%s.bed.gz | sed 's/chr//g' | bgzip -c > %s/exome_%s.bed.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("zcat %s/%s_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > %s/%s_gene_db.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("cp %s/%s_gene_names.txt.gz  %s/%s_gene_names.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    
    if annotation == True:
        
        print(
                "\n\nWARNING: The annotation step can only be performed using hg19 and hg38. Unfortunately Annovar cannot be used with grch37 and grch38. DNAscan will not perform annotation.\n\n"
            )
        
        annotation = False
    
if BED or path_gene_list:

    if path_bed:

        # splitting the analysis region into subsets of equal length to
        # distribute the work across the available threads.

        if path_gene_list:

            print(
                "\n\nWARNING: Both a bed file and a list of genes were provided. DNAscan will ignore the list of genes.\n\n"
            )

        os.system(
            "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %s > %stmp/tmp.bed"
            % (path_bed, out))

        os.system(
            "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
            % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))

        os.system("rm %stmp/tmp.bed" % (out))

        i = 0

        zero = "0"

        while i < int(num_cpu):

            if i > 9:

                zero = ""

            os.system("%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                      (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))

            #os.system("rm %stmp/%s%s" %(out,zero,str(i)))

            i += 1

    else:

        if path_gene_list:

            os.system(
                "zgrep -iwf %s %s%s_gene_names.txt.gz | awk '{print $2}' > %smatched_genes.txt"
                % (path_gene_list, path_to_db, reference, out))

            os.system(
                "zgrep -viwf %smatched_genes.txt %s  > %sunmatched_genes.txt" %
                (out, path_gene_list, out))

            if os.stat("%sunmatched_genes.txt" % (out)).st_size != 0:

                print(
                    "\n\nWARNING: some genes provided in the gene list were not found, please check which ones in %sunmatched_genes.txt "
                    % (out))

            os.system(
                "zgrep -iwf %s %s%s_gene_names.txt.gz | awk '{print $1}' > %smatched_genes_codes.txt"
                % (path_gene_list, path_to_db, reference, out))

            os.system(
                "zgrep -wf %smatched_genes_codes.txt %s%s_gene_db.txt | awk '{i=1; while (i<= int($8)) {n=split($9,a,/,/);n=split($10,b,/,/); print $2\"\t\"a[i]\"\t\"b[i]; i+=1}}' > %scustom_tmp.bed "
                % (out, path_to_db, reference, out))

            os.system(
                "%sbedtools sort -i  %scustom_tmp.bed> %scustom_sorted.bed" %
                (path_bedtools, out, out))

            os.system("%sbedtools merge -i  %scustom_sorted.bed> %scustom.bed"
                      % (path_bedtools, out, out))

            os.system("rm %scustom_sorted.bed %scustom_tmp.bed" % (out, out))

            os.system(
                "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %scustom.bed > %stmp/tmp.bed"
                % (out, out))

            os.system(
                "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
                % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))

            #os.system("rm %stmp/tmp.bed" %( out))

            i = 0

            zero = "0"

            while i < int(num_cpu):

                if i > 9:

                    zero = ""

                os.system(
                    "%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                    (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))

                # os.system("rm %stmp/%s%s" %(out,zero,str(i)))

                i += 1

        else:

            sys.exit(
                '\n\n\ERROR: the BED flag was used but neither a bed file nor a gene list was provided\n\n'
            )

#    os.system(
#        "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %
#        (path_bedtools, num_cpu, path_bed, out))
else:

    # If bed file is not provided, it generates one starting from the
    # reference genome index or the default exome bed. The pipeline uses a bed file to split the
    # analyses in several subprocesses

    if exome:

        os.system("zcat %sdb/exome_%s.bed.gz > %stmp/exome_%s.bed" %
                  (dnascan_dir, reference, out, reference))

        path_bed = "%stmp/exome_%s.bed" % (out, reference)

        if path_bed:

            # splitting the analysis region into subsets of equal length to
            # distribute the work across the available threads.

            os.system(
                "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %s > %stmp/tmp.bed"
                % (path_bed, out))

            os.system(
                "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
                % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))

            os.system("rm %stmp/tmp.bed" % (out))

            i = 0

            zero = "0"

            while i < int(num_cpu):

                if i > 9:

                    zero = ""

                os.system(
                    "%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                    (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))

                os.system("rm %stmp/%s%s" % (out, zero, str(i)))

                i += 1

    else:

        os.system("cat %s.fai | awk '{print $1\"\t0\t\"$2}' > %sreference.bed"
                  % (path_reference, out))

        path_bed = "%sreference.bed" % (out)

        os.system(
            "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'"
            % (path_bedtools, num_cpu, path_bed, out))

    BED = True

# 7. Remove duplicates command line.
# The output from HISAT2 and BWA is piped into $samblaster_cmq during the alignment
if rm_dup == "True":

    samblaster_cmq = "%ssamblaster |" % (path_samblaster)

else:

    samblaster_cmq = ""

# 8. Alignment
# Performs alignment if input sequencing data is in fastq format

if alignment:

    if format == "fastq" and "alignment.log" not in os.listdir(out + "logs"):

        # 8.1 Aligns paired end reads

        if paired == "1":

            # 8.1.1 Fast mode uses HISAT2 only to align all reads

            if mode == "fast":

                os.system(
                    "%shisat2 %s --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s  --tmpdir=%stmp -o %ssorted.bam  /dev/stdin"
                    % (path_hisat,hisat_custom_options, num_cpu, path_hisat_index, input_file,
                       input_file2, samblaster_cmq, path_samtools, num_cpu,
                       path_sambamba, num_cpu, out, out))

                bam_file = "%ssorted.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

            # 8.1.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem

            if mode == "normal" or mode == "intensive":

                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths_configs.py. Please change
                # them as needed.

                if mode == "intensive":

                    RG = True

                if RG:

                    rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)

                    rg_option_bwa = " -R '@RG\tID:%s\tLB:%s\tPL:%s\tRGPU:%s\tSM:%s' " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)

                else:

                    rg_option_hisat2 = ""

                    rg_option_bwa = ""

                os.system(
                    "%shisat2 %s %s  --no-softclip --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s  -o %ssorted.bam /dev/stdin; %ssamtools index -@ %s %ssorted.bam"
                    % (path_hisat, hisat_custom_options, rg_option_hisat2, num_cpu, path_hisat_index,
                       input_file, input_file2, samblaster_cmq, path_samtools,
                       num_cpu, path_sambamba, num_cpu, out, path_samtools,
                       num_cpu, out))

                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq"
                    % (path_samtools, num_cpu, out, out))

                os.system(
                    "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s -o %ssorted_bwa.bam  /dev/stdin ; %ssamtools index -@ %s %ssorted_bwa.bam "
                    % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                       samblaster_cmq, path_samtools, num_cpu, path_sambamba,
                       num_cpu, out, path_samtools, num_cpu, out))

                os.system("%ssamtools view -H %ssorted.bam > %sheader.txt" %
                          (path_samtools, out, out))

                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam"
                    % (path_samtools, num_cpu, out, out, out, out))

                if not debug:

                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq %sheader.txt"
                        % (out, out, out, out))

                os.system("%ssamtools index -@ %s %ssorted_merged.bam" %
                          (path_samtools, num_cpu, out))

                bam_file = "%ssorted_merged.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

        # 8.2 Aligns single end reads

        if paired == "0":

            # 8.2.1 Fast mode uses HISAT2 only to align all reads

            if mode == "fast":

                os.system(
                    "%shisat2 %s --no-spliced-alignment -p %s -x %s -U %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s -o %ssorted.bam /dev/stdin"
                    % (path_hisat, hisat_custom_options, num_cpu, path_hisat_index, input_file,
                       samblaster_cmq, path_samtools, num_cpu, path_sambamba,
                       num_cpu, out))

                bam_file = "%ssorted.bam" % (out)

                os.system("touch  %salignment.log" % (out))

            # 8.2.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem

            if mode == "normal" or mode == "intensive":

                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths_configs.py. Please change
                # them as needed.

                if mode == "intensive":

                    RG = True

                if RG:

                    rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)

                    rg_option_bwa = " -R '@RG\tID:%s\tLB:%s\tPL:%s\tRGPU:%s\tSM:%s' " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)

                else:

                    rg_option = ""

                os.system(
                    "%shisat2 %s --no-softclip --no-spliced-alignment -p %s -x %s -U %s | %s %ssamtools view -Sb -  | %ssambamba sort -t %s -o %ssorted.bam /dev/stdin; %ssamtools index -@ %s %ssorted.bam"
                    % (path_hisat, hisat_custom_options, rg_option_hisat2, num_cpu, path_hisat_index,
                       input_file, samblaster_cmq, path_samtools,
                       path_sambamba, num_cpu, out, path_samtools, num_cpu,
                       out))

                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq"
                    % (path_samtools, num_cpu, out, out))

                os.system(
                    "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq| %s %ssamtools view -Sb -  | %ssambamba sort -t %s -o %ssorted_bwa.bam /dev/stdin; %ssamtools index -@ %s %ssorted_bwa.bam "
                    % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                       samblaster_cmq, path_samtools, path_sambamba, num_cpu,
                       out, path_samtools, num_cpu, out))

                os.system("%ssamtools view -H %ssorted.bam > %sheader.txt" %
                          (path_samtools, out, out))

                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam"
                    % (path_samtools, num_cpu, out, out, out, out))

                os.system("%ssamtools index -@ %s %ssorted_merged.bam " %
                          (path_samtools, num_cpu, out))

                if not debug:

                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq "
                        % (out, out, out))

                bam_file = "%ssorted_merged.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

    else:

        if format != "fastq":

            print(
                "WARNING: Fastq format input data is requiered if you want to perform the alignment stage\n"
            )

        if "alignment.log" in os.listdir(out + "logs"):

            print(
                "WARNING: The presence of alignment.log in logs is telling you that the alignment was already peformed, please remove alignment.log if you wish to perform this stage anyway\n"
            )

            if mode == "normal" or mode == "intensive":

                bam_file = "%ssorted_merged.bam" % (out)

            if mode == "fast":

                bam_file = "%ssorted.bam" % (out)

# 9. Convert input sam file into bam

if format == "sam" and "sam2bam.log" not in os.listdir(out + "logs"):

    os.system("%ssamtools view -Sb %s  > %ssorted.bam" % (path_samtools,
                                                          input_file, out))

    bam_file = "%ssorted.bam" % (out)

    os.system("touch  %slogs/sam2bam.log" % (out))

if format == "bam":

    bam_file = "%s" % (input_file)

# 10.Variant (snv and indel) calling

if variantcalling:

    if vcf:

        print(
            "WARNING: Using input vcf as variant file. Do not provide vcf file if you wish to perform variant calling\n"
        )

        variant_results_file = vcf

    else:

        if "VC.log" in os.listdir(out + "logs"):

            print(
                "WARNING: The presence of VC.log in logs is telling you that the variant calling was already peformed, please remove VC.log if you wish to perform this stage anyway\n"
            )

            variant_results_file = "%sresults/%s_sorted.vcf.gz" % (out,
                                                                   sample_name)

        else:

            # 10.1 GATK hc indel calling (only performed in intensive mode)
            # In intesive mode, DNAscan calls snvs with Freebayes and indels
            # with GATK hc. GATK hc is used only on those positions of the
            # genome for which the alignment stage identifies one insertion or
            # deletion in at least on read.

            if mode == "intensive" and "VC_gatk.log" not in os.listdir(out +
                                                                       "logs"):

                # 10.1.1 identification of potential indel sites
                # Samtools mpileup is used to identify those positions of the genome for which the alignment stage identifies one insertion or deletion in at least on read.
                # GATK hc is used to call indels only on those positions

                counter = 1

                ps = []

                while counter < int(num_cpu) + 1:

                    command = "%ssamtools mpileup --max-depth 10000 -AQ 0 -l %stemp%s.bed %s | awk '$5 ~/[\+,\-]/ {print $1\"\t\"$2-1\"\t\"$2}' > %smpileup_positions%s.bed" % (
                        path_samtools, out, str(counter), bam_file, out,
                        str(counter))

                    proc_mpileup = subprocess.Popen(command, shell=True)

                    ps.append(proc_mpileup)

                    counter += 1

                for proc_mpileup in ps:

                    proc_mpileup.wait()

                # 10.1.2 GATK hc indel calling on selected positions (from
                # previous step 10.1.1)
                counter = 1
                os.system("touch %smpileup_positions.bed" % (out))
                while counter < int(num_cpu) + 1:

                    os.system(
                        "cat %smpileup_positions%s.bed >> %smpileup_positions.bed"
                        % (out, str(counter), out))

                    counter += 1

                os.system(
                    "%sbedtools sort -i %smpileup_positions.bed > %stemp_sorted.bed ; mv  %stemp_sorted.bed %smpileup_positions.bed"
                    % (path_bedtools, out, out, out, out))

                counter = 1

                ps = []

                while counter < int(num_cpu) + 1:

                    command = "%sjava -jar %sGenomeAnalysisTK.jar %s -R %s -T HaplotypeCaller -I %s -L %smpileup_positions%s.bed -o %sgatk_indels%s.vcf" % (
                        path_java, path_gatk, gatk_HC_custom_options, path_reference, bam_file, out,
                        str(counter), out, str(counter))

                    proc_gatk = subprocess.Popen(command, shell=True)

                    ps.append(proc_gatk)

                    counter += 1

                for proc_gatk in ps:

                    proc_gatk.wait()

                os.system(
                    "cat %sgatk_indels1.vcf | grep \"^#\" >> %sgatk_indels_merged.vcf ; for i in $(ls %s | grep gatk_indels | grep -v idx); do cat %s$i | grep -v \"^#\" >> %sgatk_indels_merged.vcf; done"
                    % (out, out, out, out, out))

                os.system(
                    "%sbedtools sort -header -i %sgatk_indels_merged.vcf > %sgatk_indels_sorted_merged.vcf"
                    % (path_bedtools, out, out))

                os.system(
                    "%svcftools  --vcf %sgatk_indels_sorted_merged.vcf --minGQ 30 --minDP 2  --recode --recode-INFO-all --out %sindels_only"
                    % (path_vcftools, out, out))

                os.system("touch  %slogs/VC_gatk.log" % (out))

            # 10.2 Freebayes snv and indel calling

            if "VC_freebayes.log" not in os.listdir(out + "logs"):

                counter = 1

                ps = []

                os.system("date")

                while counter < int(num_cpu) + 1:

                    command = path_freebayes + "freebayes " + freebayes_custom_options + " --genotype-qualities " + " -t " + out + "temp" + \
                        str(counter) + ".bed" + " -f " + path_reference + " -b " + bam_file + " > " + out + "/freebayes" + str(counter) + ".vcf"

                    proc_freebayes = subprocess.Popen(command, shell=True)

                    ps.append(proc_freebayes)

                    counter += 1

                for proc_freebayes in ps:

                    proc_freebayes.wait()

                os.system("date")

                os.system(
                    "cat %sfreebayes1.vcf | grep \"^#\" >> %smerged.vcf ; for i in $(ls %s | grep freebayes) ; do cat %s$i | grep -v \"^#\" >> %smerged.vcf; done"
                    % (out, out, out, out, out))

                os.system(
                    "%sbedtools sort -header -i %smerged.vcf > %sfreebayes.vcf"
                    % (path_bedtools, out, out))

                vcf = "freebayes.vcf"

                if mode == "intensive":

                    # If intensive mode variant calling was performed, DNAscan called snvs
                    # with Freebayes and indels with GATK hc, resulting in two vcf files.
                    # For the annotation step these two files are merged together.

                    os.system(
                        "%svcftools  --vcf %sfreebayes.vcf --minGQ 30 --minDP 2 --exclude-bed %smpileup_positions.bed  --recode --recode-INFO-all --out %sSNPs_only"
                        % (path_vcftools, out, out, out))

                    os.system("%sSNPs_only.log" % (out))

                    os.system(
                        "bgzip  %sSNPs_only.recode.vcf ; bgzip %sindels_only.recode.vcf "
                        % (out, out))

                    os.system(
                        "%stabix -p vcf %sSNPs_only.recode.vcf.gz ; %stabix -p vcf %sindels_only.recode.vcf.gz "
                        % (path_tabix, out, path_tabix, out))

                    os.system(
                        "%sjava -jar %sGenomeAnalysisTK.jar -T CombineVariants  -minimalVCF -R %s --variant %sSNPs_only.recode.vcf.gz --variant %sindels_only.recode.vcf.gz -o  %s%s.vcf --genotypemergeoption UNSORTED"
                        % (path_java, path_gatk, path_reference, out, out, out,
                           sample_name))

                    os.system(
                        "perl %svcf-sort.pl  %s%s.vcf | bgzip -c > %s%s_sorted.vcf.gz"
                        % (path_scripts, out, sample_name, out, sample_name))

                    variant_results_file = "%s%s_sorted.vcf.gz" % (out,
                                                                   sample_name)

                    os.system("%stabix -p vcf %s%s_sorted.vcf.gz" %
                              (path_tabix, out, sample_name))

                    if not debug:

                        os.system(
                            "rm %s%s.vcf* %sgatk_indels* %smpileup_positions* %sindels_only.log"
                            % (out, sample_name, out, out, out))

                else:

                    os.system(
                        "cat %sfreebayes.vcf | bgzip -c > %s%s_sorted.vcf.gz" %
                        (out, out, sample_name))

                    os.system(
                        "tabix -p vcf %s%s_sorted.vcf.gz" % (out, sample_name))

                    variant_results_file = "%s%s_sorted.vcf.gz" % (out,
                                                                   sample_name)

                if not debug:

                    os.system("rm %sfreebayes*  %stemp* %smerged.vcf" %
                              (out, out, out))

                os.system("touch  %slogs/VC_freebayes.log" % (out))

            os.system("touch  %slogs/VC.log" % (out))

# 11. Perform variant hard filtering

if filter_string:

    os.system(
        "%sbcftools filter -i \'%s\' %s | bgzip -c > %s%s_sorted_filtered.vcf.gz ; %tabix -fp vcf %s%s_sorted_filtered.vcf.gz; mv %s%s_sorted_filtered.vcf.gz %sresults " (
            path_bcftools, filter_string, variant_results_file, out,
            sample_name, path_tabix, sample_name, out, out, sample_name, out))

    variant_results_file = "%sresults/%s_sorted_filtered.vcf.gz" % (
        out, sample_name)

# 12. Performs known expansions search with ExpansionHunter
# Known expansions have to be described in a json file placed in the json
# folder (paths_configs.py).

if expansion:

    if "EH.log" in os.listdir(out):

        print(
            "WARNING: The presence of EH.log in logs is telling you that the expansion scan was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
        )

    else:

        os.system(
            "%sExpansionHunter --reads %s --reference %s  --variant-catalog %s/variant_catalog_%s.json --output-prefix %s/temp_EH"
            % (path_expansionHunter, bam_file, path_reference,
               path_expansionHunter_jsons, reference, out))

        os.system(
            "mv %s/temp_EH.vcf %s/results/%s_expansions.vcf ; bgzip %s/results/%s_expansions.vcf ; %stabix -p vcf %s/results/%s_expansions.vcf.gz"
            % (out, out, sample_name, out, sample_name, path_tabix, out,
               sample_name))

        os.system("touch  %slogs/EH.log" % (out))

        if not debug:

            os.system("rm %stemp_EH.json*  %stemp_EH.log" % (out, out))

# 13. Structural Variant calling

if SV:

    if "SV.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of SV.log in logs is telling you that structural variant calling was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
        )

    else:

        if BED:

            os.system("bgzip -c %s  > %s/temp.bed.gz" % (path_bed, out))

            os.system(
                "%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz" %
                (path_bedtools, out, out))

            os.system("%stabix -p bed %s/sorted.bed.gz" % (path_tabix, out))

            os.system("mkdir %smanta" % (out))

            os.system(
                "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta --callRegions %s/sorted.bed.gz"
                % (path_manta, bam_file, path_reference, out, out))

        else:

            os.system("mkdir %smanta" % (out))

            os.system(
                "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta"
                % (path_manta, bam_file, path_reference, out))

        os.system("%smanta/runWorkflow.py -j %s -m local" % (out, num_cpu))

        os.system(
            "mv %s/manta/results/variants/diploidSV.vcf.gz  %s/results/%s_SV.vcf.gz"
            % (out, out, sample_name))

        os.system(
            "mv %s/manta/results/variants/diploidSV.vcf.gz.tbi  %s/results/%s_SV.vcf.gz.tbi"
            % (out, out, sample_name))

        if not debug:

            os.system("rm -r %stemp.bed.gz  %ssorted.bed.gz %smanta" %
                      (out, out, out))

        os.system("touch  %slogs/SV.log" % (out))

# 14. Annotation with Annovar

if annotation:

    if "annovar.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of annovar.log in logs is telling you that annotation was already peformed, please remove annovar.log if you wish to perform this stage anyway\n"
        )

        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out,
                                                                  sample_name)

    else:

        os.system(
            "%stable_annovar.pl  --thread %s --vcfinput %s %s -buildver %s -remove -protocol %s -operation %s -nastring . --outfile %s/annovar.vcf"
            % (path_annovar, num_cpu, variant_results_file, path_annovar_db,
               reference, annovar_protocols, annovar_operations, out))
        if not debug and not alsgenescanner:
            os.system(
                "rm %sannovar.vcf.hg19_multianno.txt %sannovar.vcf.avinput" %
                (out, out))

        os.system(
            "mv %s/annovar.vcf.hg19_multianno.vcf %sresults/%s_annotated.vcf ; bgzip -f %sresults/%s_annotated.vcf ; %stabix -fp vcf %sresults/%s_annotated.vcf.gz"
            % (out, out, sample_name, out, sample_name, path_tabix, out,
               sample_name))

        os.system("mv %s %sresults/" % (variant_results_file, out))

        os.system("mv %s.tbi %sresults/" % (variant_results_file, out))

        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out,
                                                                  sample_name)

        os.system("touch  %slogs/annovar.log" % (out))

else:

    if variant_results_file:

        os.system("mv %s* %sresults/" % (variant_results_file, out))

        variant_results_file = "%sresults/%s_sorted.vcf.gz" % (out,
                                                               sample_name)

# 15. Microbes screening

if virus or bacteria or custom_microbes:

    if "microbes.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of microbes.log in logs is telling you that microbes scanning was already peformed, please remove microbes.log if you wish to perform this stage anyway\n"
        )

    else:

        # 15.1 Exctract non human reads

        os.system(
            "%ssamtools view -@ %s -hf 4 %s | %ssamtools bam2fq -s %ssingleton_reads.fastq -@ %s - > %sunaligned_reads.fastq ; cat %ssingleton_reads.fastq >> %sunaligned_reads.fastq ; gzip  %sunaligned_reads.fastq "
            % (path_samtools, num_cpu, bam_file, path_samtools, out, num_cpu,
               out, out, out, out))

        # 15.2 Identifies present non-human microbes

        if virus:

            # 15.2.1 Identifies present viruses

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -@ %s -hSb -  | %ssamtools sort -@ %s -T %stemp.file -o %soutput_virus.bam -"
                % (path_hisat, num_cpu, path_virus_index, out, path_samtools,
                   num_cpu, path_samtools, num_cpu, out, out))

            os.system(
                "%ssamtools index -@ %s %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 15.2.2 Generates virus report

            os.system(
                "awk \'{print $1}\' %svirus_stats.txt > %svirus_list.txt ; for i in $(cat %svirus_list.txt | grep C_); do printf \"$i \"; %ssamtools depth %soutput_virus.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %svirus_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            virus_coverage_stats = open('%svirus_coverage_stats.txt' % (out),
                                        'r')

            virus_coverage_stats_lines = virus_coverage_stats.readlines()

            virus_stats_file = open('%svirus_stats.txt' % (out), 'r')

            virus_stats_file_lines = virus_stats_file.readlines()

            i = 0

            virus_results = open('%sresults/virus_results.txt' % (out), 'w')

            virus_results.write("Id\tGenome_lenth\tNumber_of_reds\tCoverage\n")

            while i < len(virus_coverage_stats_lines):

                virus_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (virus_stats_file_lines[i].split('\t')[0],
                     virus_stats_file_lines[i].split('\t')[1],
                     virus_stats_file_lines[i].split('\t')[2],
                     virus_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            virus_results.close()

            if not debug:

                os.system("rm %svirus_stats.txt  %svirus_coverage_stats.txt" %
                          (out, out))

        if bacteria:

            # 15.2.3 Identifies present bacteria

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -"
                % (path_hisat, num_cpu, path_bacteria_index, out,
                   path_samtools, path_samtools, out, out))

            os.system(
                "%ssamtools index -@ %s %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 15.2.4 Generates bacteria report

            os.system(
                "awk \'{print $1}\' %sbacteria_stats.txt > %sbacteria_list.txt ; for i in $(cat %sbacteria_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_bacteria.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %sbacteria_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            bacteria_coverage_stats = open(
                '%sbacteria_coverage_stats.txt' % (out), 'r')

            bacteria_coverage_stats_lines = bacteria_coverage_stats.readlines()

            bacteria_stats_file = open('%sbacteria_stats.txt' % (out), 'r')

            bacteria_stats_file_lines = bacteria_stats_file.readlines()

            i = 0

            bacteria_results = open('%sresults/bacteria_results.txt' % (out),
                                    'w')

            bacteria_results.write(
                "Id\tGenome_lenth\tNumber_of_reds\tCoverage\n")

            while i < len(bacteria_coverage_stats_lines):

                bacteria_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (bacteria_stats_file_lines[i].split('|')[1],
                     bacteria_stats_file_lines[i].split('\t')[1],
                     bacteria_stats_file_lines[i].split('\t')[2],
                     bacteria_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            bacteria_results.close()

            if not debug:

                os.system(
                    "rm %sbacteria_stats.txt  %sbacteria_coverage_stats.txt" %
                    (out, out))

        if custom_microbes:

            # 15.2.5 Identifies present user-selected microbes

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_custom_microbes.bam -"
                % (path_hisat, num_cpu, path_custom_microbes_index, out,
                   path_samtools, path_samtools, out, out))

            os.system(
                "%ssamtools index -@ %s %soutput_custom_microbes.bam; %ssamtools idxstats %soutput_custom_microbes.bam > %scustom_microbes_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 15.2.6 Generates user-selected microbes report

            os.system(
                "awk \'{print $1}\' %scustom_microbes_stats.txt > %scustom_microbes_list.txt ; for i in $(cat %scustom_microbes_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_custom_microbes.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %scustom_microbes_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            custom_microbes_coverage_stats = open(
                '%scustom_microbes_coverage_stats.txt' % (out), 'r')

            custom_microbes_coverage_stats_lines = custom_microbes_coverage_stats.readlines(
            )

            custom_microbes_stats_file = open(
                '%scustom_microbes_stats.txt' % (out), 'r')

            custom_microbes_stats_file_lines = custom_microbes_stats_file.readlines(
            )

            i = 0

            custom_microbes_results = open(
                '%sresults/custom_microbes_results.txt' % (out), 'w')

            custom_microbes_results.write(
                "Id\tGenome_lenth\tNumber_of_reds\tCoverage\n")

            while i < len(custom_microbes_coverage_stats_lines):

                custom_microbes_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (custom_microbes_stats_file_lines[i],
                     custom_microbes_stats_file_lines[i].split('\t')[1],
                     custom_microbes_stats_file_lines[i].split('\t')[2],
                     custom_microbes_coverage_stats_lines[i].split(
                         ' ')[1].strip()))

                i += 1

            custom_microbes_results.close()

            if not debug:

                os.system(
                    "rm %scustom_microbes_stats.txt  %scustom_microbes_coverage_stats.txt"
                    % (out, out))

        os.system("touch  %slogs/microbes.log" % (out))

# 16. Alignment report generation ( samtools flagstat and stats )

if alignment_report:

    if "alignment_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of alignment_report.log in logs is telling you that the alignment report was already produced, please remove alignment_report.log if you wish to perform this stage anyway\n"
        )

    else:

        os.system("%ssamtools flagstat -@ %s %s > %sreports/%s_flagstat.txt" %
                  (path_samtools, num_cpu, bam_file, out, sample_name))

        os.system("%ssamtools stats -@ %s  %s > %sreports/%s_stats.txt" %
                  (path_samtools, num_cpu, bam_file, out, sample_name))

        os.system("touch  %slogs/alignment_report.log" % (out))

# 17. Sequencing data report generation ( fastqc )

if sequencing_report:

    if "sequencing_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of sequencing_report.log in logs is telling you that the sequence report was already produced, please remove sequencing_report.log if you wish to perform this stage anyway\n"
        )

    else:

        if path_java != "":

            java_option = "-j " + path_java + " "

        else:

            java_option = ""

        os.system("%sfastqc %s -o %sreports -f %s -t %s %s %s" %
                  (path_fastqc, java_option, out, format, num_cpu, input_file,
                   input_file2))

        os.system("touch  %slogs/sequencing_report.log" % (out))

# 18. Snv and indel calling report generation ( bcftools stats )

if calls_report:

    if "calls_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of calls_report.log in logs is telling you that the calls report was already produced, please remove calls_report.log if you wish to perform this stage anyway\n"
        )

    else:

        os.system(
            "%sbcftools stats --threads %s %s > %sreports/%s_vcfstats.txt" %
            (path_bcftools, num_cpu, variant_results_file, out, sample_name))

        os.system("touch  %slogs/calls_report.log" % (out))

# 19. Html report generation ( Multiqc )

if alignment_report or calls_report or sequencing_report:

    if "multiqc.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of multiqc.log in logs is telling you that the multiqc report was already produced, please remove multiqc.log if you wish to perform this stage anyway\n"
        )

    else:

        os.system(
            "%smultiqc -o %sreports %sreports" % (path_multiqc, out, out))

        os.system("touch  %slogs/multiqc.log" % (out))

# 20. Annotated variants report generation

if results_report:

    if "annovar.log" not in os.listdir(out + "logs") or not path_gene_list:

        print(
            "WARNING: Either the annotation was not peformed or path_gene_list was not provided in paths_configs.py, please perform annotation using the -annotation flag and specify the a gene list in paths_configs.py if you wish to perform this stage \n"
        )

    else:

        if "results_report.log" in os.listdir(out + "logs"):

            print(
                "WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n"
            )

        else:

            os.system("zcat %s > %stemp.vcf" % (variant_results_file, out))

            #vcf = open('%stemp.vcf' % (out), 'r')

            #vcf_lines = vcf.readlines()

            gene_list_file = open(path_gene_list)

            gene_list_lines = gene_list_file.readlines()

            gene_list = gene_list_lines

            out_file_all = open(
                '%sreports/%s_all_variants.txt' % (out, sample_name), 'w')

            counter = 0

            for i in gene_list:

                with open('%stemp.vcf' % (out)) as vcf:

                    for j in vcf:

                        check1 = re.search(
                            r'(^chr)|(^[0-9,X,Y,M]+\t)', j, flags=0)

                        check = re.search(
                            "=%s;" % (i.strip().upper()), j, flags=0)

                        if check and check1:

                            infos = j.split('ANNOVAR_DATE')[1][12:].split(
                                'ALLELE_END')[0].replace(";", "\t")

                            if counter == 0:

                                replaced_1 = re.sub(
                                    '=[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+', '',
                                    infos)

                                out_file_all.write(
                                    'CHR\tPosition\tRef\tAlt\tGenotype\t%s\n' %
                                    (replaced_1))

                                counter = 1

                            replaced = re.sub('[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+=',
                                              '', infos)

                            out_file_all.write(
                                '%s\t%s\t%s\t%s\t%s\t%s\n' %
                                (j.split('\t')[0], j.split('\t')[1],
                                 j.split('\t')[3], j.split('\t')[4],
                                 j.split('\t')[-1].split(':')[0], replaced))

            out_file_all.close()

            os.system("touch  %slogs/results_report.log" % (out))

# 21. Starting iobio services

if iobio:
    if "iobio.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of iobio.log in logs is telling you that the iobio services were already started, please remove iobio.log if you wish to start them again\n"
        )

    else:

        os.system("cd %s ; nohup python3 -m http.server %s  >/dev/null 2>&1 &"
                  % (path_iobio, port_num))

        os.system("touch  %slogs/iobio.log" % (out))

        print(
            "\n\nIobio serces have been started at http://localhost:%s\n\nCopy and paste http://localhost:%s to select the service (vcf, bam, gene) and upload your data into the selected service\n\nIf you want to explore your variant calling results please copy and paste the following URL into your browser and upload the vcf file (../%sresults/%s_sorted.vcf.gz):\n\n"
            % (port_num, port_num, out, sample_name),
            end='',
            flush=True)

        if "annovar.log" in os.listdir(out + "logs"):

            print(
                "http://localhost:%s/gene.iobio/?species=Human&rel0=proband&rel1=mother&rel2=father&genes="
                % (port_num),
                end='',
                flush=True)

            a = {}

            os.system("zcat %s > %stemp.vcf" % (variant_results_file, out))

            file = open('%stemp.vcf' % (out), 'r')

            file_lines = file.readlines()

            for i in file_lines:

                check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

                if check:

                    a[i.split('Gene.refGene=')[1].split(';')[0]] = []

            for i in a.keys():

                print('%s,' % (i.split(',')[0]), end='', flush=True)
            print(' \n\n')

        else:

            print("http://localhost:%s/gene.iobio/" % (port_num))

if alsgenescanner:

    os.system(
        "python3 %s/alsgenescanner.py %s/annovar.vcf.hg19_multianno.txt %s/results/%s_alsgenescanner_all.txt"
        % (path_scripts, out, out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_alsod.txt; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_alsod.txt >> %s/results/%s_alsgenescanner_alsod.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_clinvar.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_clinvar.txt >> %s/results/%s_alsgenescanner_clinvar.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_manual_review.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_manual_review.txt >> %s/results/%s_alsgenescanner_manual_review.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_all_ranked.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -i pathog | sed 's/ /_/g'| sort -k10nr >> %s/results/%s_alsgenescanner_all_ranked.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep '^chr'  | grep -iv pathog | sed 's/ /_/g'| sort -k10nr >> %s/results/%s_alsgenescanner_all_ranked.txt  "
        % (out, sample_name, out, sample_name, out, sample_name, out, sample_name, 
           out, sample_name, out, sample_name))
    
