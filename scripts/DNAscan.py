################################################################
# Program: DNAscan
# Version 1.0
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################


################################################################
# Script structure:
# 1. Import python modules
# 2. Define paths viriables
# 3. Define options from command line
# 4. Parse options from command line
# 5. Create working dir tree
# 6. Bed splitting
# 7. Alignment
#   7.1 Aligns paired end reads
#       7.1.1 Fast mode alignment
#       7.1.2 Normal and intesive mode alignment
#   7.2 Aligns single end reads
#       7.2.1 Fast mode alignment
#       7.2.2 Normal and intesive mode alignment
# 8. Remove duplicates
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
# 18. Snv and indel calling report generation ( rtg vcfstats )
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
import paths
import subprocess

from argparse import RawTextHelpFormatter


# 2. Define paths variables from paths.py

path_iobio = paths.path_iobio

port_num = paths.port_num

path_vcftools = paths.path_vcftools

path_gatk = paths.path_gatk

path_java = paths.path_java

path_hisat = paths.path_hisat

path_hisat_index = paths.path_hisat_index

path_samtools = paths.path_samtools

path_freebayes = paths.path_freebayes

path_annovar = paths.path_annovar

path_annovar_db = paths.path_annovar_db

num_cpu = paths.num_cpu

path_bed = paths.path_bed

path_expansionHunter = paths.path_expansionHunter

path_reference = paths.path_reference

path_expansionHunter_jsons = paths.path_expansionHunter_jsons

path_manta = paths.path_manta

path_bedtools = paths.path_bedtools

path_tabix = paths.path_tabix

path_virus_index = paths.path_virus_index

path_bacteria_index = paths.path_virus_index

path_custom_microbes_index = paths.path_custom_microbes_index

path_bwa = paths.path_bwa

path_bwa_index = paths.path_bwa_index

path_rtg = paths.path_rtg

path_bcftools = paths.path_bcftools

RG_ID = paths.RG_ID

RG_LB = paths.RG_LB

RG_PL = paths.RG_PL

RG_PU = paths.RG_PU

RG_SM = paths.RG_SM

path_scripts = paths.path_scripts


# 3. Define options variables from command line

parser = argparse.ArgumentParser(
    prog='python DNAscan.py',
    usage='%(prog)s [options] -format "string" -reference "string" -in "string"',
    description='############DNAscan help Message############ \n\nDNAscan uses the file paths.py to locate the needed tools and files. Please make sure your paths.py file is properly filled \n\nUsage example: \n\npython rDNAscan.py -format fastq -out /home/user/test/ -in sample_sorted_1.fq.gz -reference hg19 -alignment -variantcalling -results_report\n\nPlease check the following list of optional and required options\n\n################################################',
    formatter_class=RawTextHelpFormatter)

parser.add_argument(
    '-RG',
    action="store",
    dest="RG",
    default=False,
    help='if this flag is set the alignment stage will add the provided in paths.py read group (Default = "False")')


parser.add_argument(
    '-mode',
    required=True,
    action="store",
    dest="mode",
    default="fast",
    help='options are fast, normal, intensive [string] (default = "fast")')

parser.add_argument(
    '-filter_string',
    action="store",
    dest="filter_string",
    default="",
    help='bcftools filter string, eg GQ>20 & DP>10 (Default = "")')

parser.add_argument(
    '-paired',
    action="store",
    dest="paired",
    default="1",
    help='options are 1 for paired end reads and 0 for single end reads (Default = "1")')

parser.add_argument(
    '-vcf',
    action="store",
    dest="vcf",
    help='complementary vcf file')

parser.add_argument(
    '-in2',
    action="store",
    dest="input_file2",
    help='input file 2, for paired end reads only (usually fastq file)')

parser.add_argument(
    '-iobio',
    action="store_true",
    dest="iobio",
    help='if this flag is set iobio services will be started at the end of the analysis (Default = "False")',
    default=False)

parser.add_argument(
    '-alignment',
    action="store_true",
    dest="alignment",
    help='if this flag is set the alignment stage will be performed (Default = "False")',
    default=False)

parser.add_argument(
    '-expansion',
    action="store_true",
    dest="expansion",
    help='if this flag is set DNAscan will look for the expansions described in the jason folder described in paths.py  (Default = "False") ',
    default=False)

parser.add_argument(
    '-SV',
    action="store_true",
    dest="SV",
    help='if this flag is set the structural variant calling stage will be performed (Default = "False") ',
    default=False)

parser.add_argument(
    '-BED',
    action="store_true",
    dest="BED",
    help='restrict the analysis to the regions in the bed file (Default = "False") ',
    default=False)

parser.add_argument(
    '-virus',
    action="store_true",
    dest="virus",
    help='if this flag is set DNAscan will perform viral scanning (Default = "False")  ',
    default=False)

parser.add_argument(
    '-bacteria',
    action="store_true",
    dest="bacteria",
    help='if this flag is set DNAscan will perform bacteria scanning (Default = "False") ',
    default=False)

parser.add_argument(
    '-custom_microbes',
    action="store_true",
    dest="custom_microbes",
    help='if this flag is set DNAscan will perform a customized microbe scanning according to the provided microbe data base in paths.py (Default = "False")  ',
    default=False)

parser.add_argument(
    '-variantcalling',
    action="store_true",
    dest="variantcalling",
    help='if this flag is set DNAscan will perform snv and indel calling (Default = "False")  ',
    default=False)

parser.add_argument(
    '-annotation',
    action="store_true",
    dest="annotation",
    help='if this flag is set DNAscan will annotate the found variants (Default = "False")  ',
    default=False)

parser.add_argument(
    '-results_report',
    action="store_true",
    dest="results_report",
    help='if this flag is set DNAscan generate a results report (Default = "False")  ',
    default=False)

parser.add_argument(
    '-alignment_report',
    action="store_true",
    dest="alignment_report",
    help='if this flag is set DNAscan generate an alignment report (Default = "False") ',
    default=False)

parser.add_argument(
    '-sequencing_report',
    action="store_true",
    dest="sequencing_report",
    help='if this flag is set DNAscan generate a report describing the input sequencing data (Default = "False") ',
    default=False)

parser.add_argument(
    '-calls_report',
    action="store_true",
    dest="calls_report",
    help='if this flag is set DNAscan generate a report describing the found snvs and indels (Default = "False")',
    default=False)

parser.add_argument(
    '-sample_name',
    action="store",
    dest="sample_name",
    default="sample",
    help='specify sample name [string] (default = "sample")')

parser.add_argument(
    '-rm_dup',
    action="store_true",
    dest="rm_dup",
    help='if this flag is set DNAscan will remove duplicates from the sample bam file (Default = "False")',
    default=False)

parser.add_argument(
    '-debug',
    action="store_true",
    dest="debug",
    help='if this flag is set DNAscan will not delete intermediete and temporary files (Default = "False")',
    default=False)

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument(
    '-out',
    required=True,
    action="store",
    dest="out",
    help='path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/')

requiredNamed.add_argument(
    '-format',
    required=True,
    action="store",
    dest="format",
    default="fastq",
    help='options are bam, sam, fastq, vcf [string] ')

requiredNamed.add_argument(
    '-reference',
    required=True,
    action="store",
    dest="reference",
    help='options are hg19, hg38 [string]')

requiredNamed.add_argument(
    '-in',
    required=True,
    action="store",
    dest="input_file",
    help='input file [string]')


# 4. Parse command line options

args = parser.parse_args()

mode = args.mode

format = args.format

paired = args.paired

reference = args.reference

input_file = args.input_file

input_file2 = args.input_file2

expansion = args.expansion

out = args.out

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

# 5. Create working dir tree

os.system(
    "mkdir  %slogs ; mkdir  %sreports ; mkdir  %sresults" %
    (out, out, out))


# 6. Bed splitting: splitting the analysis region into subsents of equal length to distribute the work across the available threads.
# To do this DNAscan uses a bed file.
# If bed file is not provided, it generates one starting from the
# reference genome index. The pipeline uses a bed file to split the
# analyses in several subprocesses


if BED:

    # splitting the analysis region into subsents of equal length to
    # distribute the work across the available threads.

    os.system(
        "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %
        (path_bedtools, num_cpu, path_bed, out))

else:

    # If bed file is not provided, it generates one starting from the
    # reference genome index. The pipeline uses a bed file to split the
    # analyses in several subprocesses

    os.system(
        "cat %s.fai | awk '{print $1\"\t0\t\"$2}' > %sreference.bed" %
        (path_reference, out))

    path_bed = "%sreference.bed" % (out)

    os.system(
        "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %
        (path_bedtools, num_cpu, path_bed, out))

    BED = True

# 7. Alignment
# Performs alignment if input sequencing data is in fastq format

if alignment:

    if format == "fastq" and "alignment.log" not in os.listdir(out + "logs"):

        # 7.1 Aligns paired end reads

        if paired == "1":

            # 7.1.1 Fast mode uses HISAT2 only to align all reads

            if mode == "fast":

                os.system(
                    "%shisat2 --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %ssamtools view -Sb -  | %ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %
                    (path_hisat,
                     num_cpu,
                     path_hisat_index,
                     input_file,
                     input_file2,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                bam_file = "%ssorted.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

            # 7.1.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem

            if mode == "normal" or mode == "intensive":

                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths.py. Please change
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
                    "%shisat2 %s  --no-softclip --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted.bam ; %ssamtools index %ssorted.bam" %
                    (path_hisat,
                     rg_option_hisat2,
                     num_cpu,
                     path_hisat_index,
                     input_file,
                     input_file2,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq" %
                    (path_samtools, num_cpu, out, out))

                os.system(
                    "%sbwa mem %s -t %s %s %sunaligned_reads.fq | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted_bwa.bam ; %ssamtools index %ssorted_bwa.bam " %
                    (path_bwa,
                     rg_option_bwa,
                     num_cpu,
                     path_bwa_index,
                     out,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                os.system(
                    "%ssamtools view -H %ssorted.bam > %sheader.txt" %
                    (path_samtools, out, out))

                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam" %
                    (path_samtools, num_cpu, out, out, out, out))

                if not debug:

                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq %sheader.txt" %
                        (out, out, out, out))

                os.system(
                    "%ssamtools index %ssorted_merged.bam" %
                    (path_samtools, out))

                bam_file = "%ssorted_merged.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

        # 7.2 Aligns single end reads

        if paired == "0":

            # 7.2.1 Fast mode uses HISAT2 only to align all reads

            if mode == "fast":

                os.system(
                    "%shisat2 --no-spliced-alignment--remove-chrname -p %s -x %s -U %s | %ssamtools view -Sb -  | %ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %
                    (path_hisat,
                     num_cpu,
                     path_hisat_index,
                     input_file,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                bam_file = "%ssorted.bam" % (out)

                os.system("touch  %salignment.log" % (out))

            # 7.2.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem

            if mode == "normal" or mode == "intensive":

                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths.py. Please change
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
                    "%shisat2  --no-softclip --no-spliced-alignment -p %s -x %s -U %s | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted.bam ; %ssamtools index %ssorted.bam" %
                    (path_hisat,
                     rg_option_hisat2,
                     num_cpu,
                     path_hisat_index,
                     input_file,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq" %
                    (path_samtools, num_cpu, out, out))

                os.system(
                    "%sbwa mem %s -t %s %s %sunaligned_reads.fq| %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted_bwa.bam ; %ssamtools index %ssorted_bwa.bam " %
                    (path_bwa,
                     rg_option_bwa,
                     num_cpu,
                     path_bwa_index,
                     out,
                     path_samtools,
                     path_samtools,
                     out,
                     out,
                     path_samtools,
                     out))

                os.system(
                    "%ssamtools view -H %ssorted.bam > %sheader.txt" %
                    (path_samtools, out, out))

                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam" %
                    (path_samtools, num_cpu, out, out, out, out))

                os.system(
                    "%ssamtools index %ssorted_merged.bam " %
                    (path_samtools, out))

                if not debug:

                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq " %
                        (out, out, out))

                bam_file = "%ssorted_merged.bam" % (out)

                os.system("touch  %slogs/alignment.log" % (out))

    else:

        if format != "fastq":

            print(
                "WARNING: Fastq format input data is requiered if you want to perform the alignment stage\n")

        if "alignment.log" in os.listdir(out + "logs"):

            print("WARNING: The presence of alignment.log in logs is telling you that the alignment was already peformed, please remove alignment.log if you wish to perform this stage anyway\n")

            if mode == "normal" or mode == "intensive":

                bam_file = "%ssorted_merged.bam" % (out)

            if mode == "fast":

                bam_file = "%ssorted.bam" % (out)

# 8. Remove duplicates

if rm_dup:

    if "rm_dup.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of rm_dup.log in logs is telling you that the duplicate removal was already peformed, please remove rm_dup.log if you wish to perform this stage anyway\n")

    else:

        if paired == 0:

            single_end_option = "-S"

        else:

            single_end = ""

        os.system("%ssamtools rmdup %s %s %ssorted_merged_rmdup.bam" %
                  (path_samtools, single_end_option, bam_file, out))

        bam_file = "%ssorted_merged_rmdup.bam" % (out)

        if not debug:

            os.system("rm %ssorted_merged.bam" % (out))

        s.system("touch  %slogs/rm_dup.log" % (out))

# 9. Convert input sam file into bam

if format == "sam" and "sam2bam.log" not in os.listdir(out + "logs"):

    os.system(
        "%ssamtools view -Sb %s  > %ssorted.bam" %
        (path_samtools, input_file, out))

    bam_file = "%ssorted.bam" % (out)

    os.system("touch  %slogs/sam2bam.log" % (out))


# 10.Variant (snv and indel) calling

if variantcalling:

    if vcf:

        print("WARNING: Using input vcf as variant file. Do not provide vcf file if you wish to perform variant calling\n")
        
        variant_results_file = vcf

    else:

        if "VC.log" in os.listdir(out + "logs"):

            print("WARNING: The presence of VC.log in logs is telling you that the variant calling was already peformed, please remove VC.log if you wish to perform this stage anyway\n")

            variant_results_file = "%sresults/%s_sorted.vcf.gz" % (out, sample_name)

        else:

            # 10.1 GATK hc indel calling (only performed in intensive mode)
            # In intesive mode, DNAscan calls snvs with Freebayes and indels
            # with GATK hc. GATK hc is used only on those positions of the
            # genome for which the alignment stage identifies one insertion or
            # deletion in at least on read.

            if mode == "intensive" and "VC_gatk.log" not in os.listdir(
                    out + "logs"):

                # 10.1.1 identification of potential indel sites
                # Samtools mpileup is used to identify those positions of the genome for which the alignment stage identifies one insertion or deletion in at least on read.
                # GATK hc is used to call indels only on those positions

                counter = 1

                ps = []

                while counter < int(num_cpu) + 1:

                    command = "%ssamtools mpileup --max-depth 10000 -AQ 0 -l %stemp%s.bed %s | awk '$5 ~/[\+,\-]/ {print $1\"\t\"$2-1\"\t\"$2}' > %smpileup_positions%s.bed" % (
                        path_samtools, out, str(counter), bam_file, out, str(counter))

                    proc_mpileup = subprocess.Popen(command, shell=True)

                    ps.append(proc_mpileup)

                    counter += 1

                for proc_mpileup in ps:

                    proc_mpileup.wait()

                # 10.1.2 GATK hc indel calling on selected positions (from
                # previous step 10.1.1)

                counter = 1

                ps = []

                while counter < int(num_cpu) + 1:

                    command = "%sjava -jar %sGenomeAnalysisTK.jar -R %s -T HaplotypeCaller -I %s -L %smpileup_positions%s.bed -o %sgatk_indels%s.vcf" % (
                        path_java, path_gatk, path_reference, bam_file, out, str(counter), out, str(counter))

                    proc_gatk = subprocess.Popen(command, shell=True)

                    ps.append(proc_gatk)

                    counter += 1

                for proc_gatk in ps:

                    proc_gatk.wait()

                os.system(
                    "cat %sgatk_indels1.vcf | grep \"^#\" >> %sgatk_indels_merged.vcf ; for i in $(ls %s | grep gatk_indels | grep -v idx); do cat %s$i | grep -v \"^#\" >> %sgatk_indels_merged.vcf; done" %
                    (out,
                     out,
                     out,
                     out,
                     out))

                os.system(
                    "%sbedtools sort -header -i %sgatk_indels_merged.vcf > %sgatk_indels_sorted_merged.vcf" %
                    (path_bedtools, out, out))

                os.system(
                    "%svcftools  --vcf %sgatk_indels_sorted_merged.vcf --minGQ 20 --minDP 2 --keep-only-indels --recode --recode-INFO-all --out %sindels_only" %
                    (path_vcftools, out, out))

                if not debug:

                    os.system(
                        "rm %sgatk_indels* %smpileup_positions* %sindels_only.log" %
                        (out, out, out))

                os.system("touch  %slogs/VC_gatk.log" % (out))

            # 10.2 Freebayes snv and indel calling

            if "VC_freebayes.log" not in os.listdir(out + "logs"):

                counter = 1

                ps = []

                os.system("date")

                while counter < int(num_cpu) + 1:

                    command = path_freebayes + "freebayes" + " -t " + out + "temp" + \
                        str(counter) + ".bed" + " -f " + path_reference + " -b " + bam_file + " > " + out + "/freebayes" + str(counter) + ".vcf"

                    proc_freebayes = subprocess.Popen(command, shell=True)

                    ps.append(proc_freebayes)

                    counter += 1

                for proc_freebayes in ps:

                    proc_freebayes.wait()

                os.system("date")

                os.system(
                    "cat %sfreebayes1.vcf | grep \"^#\" >> %smerged.vcf ; for i in $(ls %s | grep freebayes) ; do cat %s$i | grep -v \"^#\" >> %smerged.vcf; done" %
                    (out, out, out, out, out))

                os.system(
                    "%sbedtools sort -header -i %smerged.vcf > %sfreebayes.vcf" %
                    (path_bedtools, out, out))

                vcf = "freebayes.vcf"

                if mode == "intensive":


                    # If intensive mode variant calling was performed, DNAscan called snvs
                    # with Freebayes and indels with GATK hc, resulting in two vcf files.
                    # For the annotation step these two files are merged together.
                    
                    
                    os.system(
                        "%svcftools  --vcf %sfreebayes.vcf --minGQ 20 --minDP 2 --remove-indels --recode --recode-INFO-all --out %sSNPs_only" %
                        (path_vcftools, out, out))

                    os.system(
                        "bgzip  %sSNPs_only.recode.vcf ; bgzip %sindels_only.recode.vcf " %
                        (out, out))

                    os.system(
                        "%stabix -p vcf %sSNPs_only.recode.vcf.gz ; %stabix -p vcf %sindels_only.recode.vcf.gz " %
                        (path_tabix, out, path_tabix, out))
                        
                    os.system(
                        "%sjava -jar %sGenomeAnalysisTK.jar -T CombineVariants -R %s --variant %sSNPs_only.recode.vcf.gz --variant %sindels_only.recode.vcf.gz -o  %s%s.vcf --genotypemergeoption UNSORTED" %
                        (path_java,
                         path_gatk,
                         path_reference,
                         out,
                         out,
                         out,                      
                         sample_name))

                    os.system(
                        "perl %svcf-sort.pl  %s%s.vcf | bgzip -c > %s%s_sorted.vcf.gz" %
                        (path_scripts, out, sample_name, out, sample_name))
                    
                    variant_results_file = "%s%s_sorted.vcf.gz" %(out, sample_name)

                    

                    os.system(
                            "%stabix -p vcf %s%s_sorted.vcf.gz" %
                            (path_tabix, out, sample_name))
                    
                    if not debug:

                            os.system("rm %s%s.vcf* " % (out, sample_name))

                else:

                    os.system(
                        "cat %sfreebayes.vcf | bgzip -c > %s%s_sorted.vcf.gz" %
                        (out, out, sample_name))

                    os.system(
                        "tabix -p vcf %s%s_sorted.vcf.gz" %
                        (out, sample_name))


                    variant_results_file = "%s%s_sorted.vcf.gz" %(out, sample_name)

                if not debug:

                    os.system(
                        "rm %sfreebayes* %sSNPs_only.log %stemp* %smerged.vcf" %
                        (out, out, out, out))

                os.system("touch  %slogs/VC_freebayes.log" % (out))

            os.system("touch  %slogs/VC.log" % (out))


# 11. Perform variant hard filtering

if filter_string:
    
    
    os.system(
        "%sbcftools filter -i \'%s\' %s | bgzip -c > %s%s_sorted_filtered.vcf.gz ; %tabix -fp vcf %s%s_sorted_filtered.vcf.gz" (
            path_bcftools,
            filter_string,
            variant_results_file,
            out,
            sample_name,
            path_tabix,
            sample_name,
            out))
            
    variant_results_file = "%s%s_sorted_filtered.vcf.gz" %(out, sample_name)
''''

    if mode == "intensive":

        os.system(
            "%sbcftools filter -i \'%s\' %sSNPs_only.recode.vcf.gz | bgzip -c > %stemp1.vcf.gz ; mv %stemp1.vcf.gz %sSNPs_only.recode.vcf.gz ; %tabix -fp vcf %sSNPs_only.recode.vcf.gz" (
                path_bcftools,
                filter_string,
                out,
                out,
                out,
                out,
                path_tabix,
                out))

        os.system(
            "%sbcftools filter -i \'%s\' %sindels_only.recode.vcf.gz | bgzip -c > %stemp2.vcf.gz ; mv %stemp2.vcf.gz %sindels_only.recode.vcf.gz ; %tabix -fp vcf %sindels_only.recode.vcf.gz" (
                path_bcftools,
                filter_string,
                out,
                out,
                out,
                out,
                path_tabix,
                out))

    else:

        os.system(
            "%sbcftools filter -i \'%s\' %s%s_sorted.vcf.gz| bgzip -c > %stemp1.vcf.gz ; mv %stemp1.vcf.gz %s%s_sorted.vcf.gz ; %tabix -fp vcf %s%s_sorted.vcf.gz" (
                path_bcftools,
                filter_string,
                out,
                sample_name,
                out,
                out,
                sample_name,
                out,
                path_tabix,
                out,
                sample_name))

'''
# 12. Performs known expansions search with ExpansionHunter
# Known expansions have to be described in a json file placed in the json
# folder (paths.py).

if expansion:

    if "EH.log" in os.listdir(out):

        print("WARNING: The presence of EH.log in logs is telling you that the expansion scan was already peformed, please remove SV.log if you wish to perform this stage anyway\n")

    else:

        os.system(
            "%sExpansionHunter --bam %s --ref-fasta %s --repeat-specs %s --vcf %s/EH_vcf --json %s/test.json --log %s/test_log" %
            (path_expansionHunter, bam_file, path_reference, path_expansionHunter_jsons, out, out, out))

        os.system(
            "mv %s/EH_vcf %s/results/%s_expansions.vcf ; bgzip %s/results/%s_expansions.vcf ; %stabix -p vcf %s/results/%s_expansions.vcf.gz" %
            (out, out, sample_name, out, sample_name, path_tabix, out, sample_name))

        os.system("touch  %slogs/EH.log" % (out))


# 13. Structural Variant calling

if SV:

    if "SV.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of SV.log in logs is telling you that structural variant calling was already peformed, please remove SV.log if you wish to perform this stage anyway\n")

    else:

        if BED:

            os.system("bgzip -c %s  > %s/temp.bed.gz" % (path_bed, out))

            os.system(
                "%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz" %
                (path_bedtools, out, out))

            os.system("%stabix -p bed %s/sorted.bed.gz" % (path_tabix, out))

            os.system("mkdir %smanta" % (out))

            os.system(
                "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta --callRegions %s/sorted.bed.gz" %
                (path_manta, bam_file, path_reference, out, out))

        else:

            os.system("mkdir %smanta" % (out))

            os.system(
                "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta" %
                (path_manta, bam_file, path_reference, out))

        os.system("%smanta/runWorkflow.py -j %s -m local" % (out, num_cpu))

        os.system(
            "mv %s/manta/results/variants/diploidSV.vcf.gz  %s/results/%s_SV.vcf.gz" %
            (out, out, sample_name))

        os.system(
            "mv %s/manta/results/variants/diploidSV.vcf.gz.tbi  %s/results/%s_SV.vcf.gz.tbi" %
            (out, out, sample_name))

        if not debug:

            os.system(
                "rm -r %stemp.bed.gz  %ssorted.bed.gz %smanta" %
                (out, out, out))

        os.system("touch  %slogs/SV.log" % (out))

# 14. Annotation with Annovar

if annotation:

    if "annovar.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of annovar.log in logs is telling you that annotation was already peformed, please remove annovar.log if you wish to perform this stage anyway\n")
        
        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out, sample_name)
   
    else:
    
        os.system(
            "%stable_annovar.pl  --thread %s --vcfinput %s %s -buildver %s -remove -protocol refGene,dbnsfp30a,clinvar_20170130,avsnp147,cadd -operation g,f,f,f,f -nastring . --outfile %s/annovar.vcf" %
            (path_annovar,
             num_cpu,
             variant_results_file,
             path_annovar_db,
             reference,
             out))

        os.system(
            "rm %sannovar.vcf.hg19_multianno.txt %sannovar.vcf.avinput" %
            (out, out))
        
        os.system(
            "mv %s/annovar.vcf.hg19_multianno.vcf %sresults/%s_annotated.vcf ; bgzip %sresults/%s_annotated.vcf ; %stabix -fp vcf %sresults/%s_annotated.vcf.gz" %
            (out,
             out,
             sample_name,
             out,
             sample_name,
             path_tabix,
             out,
             sample_name))
        
        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out, sample_name)
        
        os.system("touch  %slogs/annovar.log" % (out))
        
else :
    
    os.system("mv %s* %sresults/" %(variant_results_file, out))
    
    variant_results_file = "%sresults/%s_sorted.vcf.gz" % (out, sample_name)
    


# 15. Microbes screening

if virus or bacteria or custom_microbes:

    if "microbes.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of microbes.log in logs is telling you that microbes scanning was already peformed, please remove microbes.log if you wish to perform this stage anyway\n")

    else:

        # 15.1 Exctract non human reads

        os.system(
            "%ssamtools view -hf 4 %s | %ssamtools bam2fq - | gzip -c - > %suniligned_reads.fastq.gz" %
            (path_samtools, bam_file, path_samtools, out))

        # 15.2 Identifies present non-human microbes

        if virus:

            # 15.2.1 Identifies present viruses

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_virus.bam -" %
                (path_hisat, num_cpu, path_virus_index, out, path_samtools, path_samtools, out, out))

            os.system(
                "%ssamtools index %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt" %
                (path_samtools, out, path_samtools, out, out))

            # 15.2.2 Generates virus report

            os.system(
                "awk \'{print $1}\' %svirus_stats.txt > %svirus_list.txt ; for i in $(cat %svirus_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_virus.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %svirus_coverage_stats.txt" %
                (out,
                 out,
                 out,
                 path_samtools,
                 out,
                 out))

            virus_coverage_stats = open(
                '%svirus_coverage_stats.txt' %
                (out), 'r')

            virus_coverage_stats_lines = virus_coverage_stats.readlines()

            virus_stats_file = open('%svirus_stats.txt' % (out), 'r')

            virus_stats_file_lines = virus_stats_file.readlines()

            i = 0

            virus_results = open('%sresults/virus_results.txt' % (out), 'w')

            virus_results.write("Id\tGenome_lenth\tNumber_of_reds\tCoverage\n")

            while i < len(virus_coverage_stats_lines):

                virus_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (virus_stats_file_lines[i].split('|')[1],
                     virus_stats_file_lines[i].split('\t')[1],
                     virus_stats_file_lines[i].split('\t')[2],
                     virus_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            virus_results.close()

            if not debug:

                os.system(
                    "rm %svirus_stats.txt  %svirus_coverage_stats.txt" %
                    (out, out))

        if bacteria:

            # 15.2.3 Identifies present bacteria

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -" %
                (path_hisat, num_cpu, path_bacteria_index, out, path_samtools, path_samtools, out, out))

            os.system(
                "%ssamtools index %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt" %
                (path_samtools, out, path_samtools, out, out))

            # 15.2.4 Generates bacteria report

            os.system(
                "awk \'{print $1}\' %sbacteria_stats.txt > %sbacteria_list.txt ; for i in $(cat %sbacteria_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_bacteria.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %sbacteria_coverage_stats.txt" %
                (out,
                 out,
                 out,
                 path_samtools,
                 out,
                 out))

            bacteria_coverage_stats = open(
                '%sbacteria_coverage_stats.txt' %
                (out), 'r')

            bacteria_coverage_stats_lines = bacteria_coverage_stats.readlines()

            bacteria_stats_file = open('%sbacteria_stats.txt' % (out), 'r')

            bacteria_stats_file_lines = bacteria_stats_file.readlines()

            i = 0

            bacteria_results = open(
                '%sresults/bacteria_results.txt' %
                (out), 'w')

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
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_custom_microbes.bam -" %
                (path_hisat,
                 num_cpu,
                 path_custom_microbes_index,
                 out,
                 path_samtools,
                 path_samtools,
                 out,
                 out))

            os.system(
                "%ssamtools index %soutput_custom_microbes.bam; %ssamtools idxstats %soutput_custom_microbes.bam > %scustom_microbes_stats.txt" %
                (path_samtools, out, path_samtools, out, out))

            # 15.2.6 Generates user-selected microbes report

            os.system(
                "awk \'{print $1}\' %scustom_microbes_stats.txt > %scustom_microbes_list.txt ; for i in $(cat %scustom_microbes_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_custom_microbes.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %scustom_microbes_coverage_stats.txt" %
                (out,
                 out,
                 out,
                 path_samtools,
                 out,
                 out))

            custom_microbes_coverage_stats = open(
                '%scustom_microbes_coverage_stats.txt' %
                (out), 'r')

            custom_microbes_coverage_stats_lines = custom_microbes_coverage_stats.readlines()

            custom_microbes_stats_file = open(
                '%scustom_microbes_stats.txt' %
                (out), 'r')

            custom_microbes_stats_file_lines = custom_microbes_stats_file.readlines()

            i = 0

            custom_microbes_results = open(
                '%sresults/custom_microbes_results.txt' %
                (out), 'w')

            custom_microbes_results.write(
                "Id\tGenome_lenth\tNumber_of_reds\tCoverage\n")

            while i < len(custom_microbes_coverage_stats_lines):

                custom_microbes_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (custom_microbes_stats_file_lines[i],
                     custom_microbes_stats_file_lines[i].split('\t')[1],
                     custom_microbes_stats_file_lines[i].split('\t')[2],
                     custom_microbes_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            custom_microbes_results.close()

            if not debug:

                os.system(
                    "rm %scustom_microbes_stats.txt  %scustom_microbes_coverage_stats.txt" %
                    (out, out))

        os.system("touch  %slogs/microbes.log" % (out))


# 16. Alignment report generation ( samtools flagstat and stats )


if alignment_report:

    if "alignment_report.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of alignment_report.log in logs is telling you that the alignment report was already produced, please remove alignment_report.log if you wish to perform this stage anyway\n")

    else:

        os.system(
            "%ssamtools flagstat -@ %s %s > %sreport/%s_flagstat.log" %
            (path_samtools, num_cpu, bam_file, out, sample_name))

        os.system("%ssamtools stats -@ %s %s > %sreport/%s_stats.log" %
                  (path_samtools, num_cpu, bam_file, out, sample_name))

        os.system("touch  %slogs/alignment_report.log" % (out))

# 17. Sequencing data report generation ( fastqc )

if sequencing_report:

    if "sequencing_report.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of sequencing_report.log in logs is telling you that the sequence report was already produced, please remove sequencing_report.log if you wish to perform this stage anyway\n")

    else:

        if path_java != "":

            java_option = "-j " + path_java + " "

        os.system("%sfastqc %s -o %sreports -f %s -t %s %s" %
                  (path_fastqc, java_option, out, format, num_cpu, nput_file))

        os.system("touch  %slogs/sequencing_report.log" % (out))


# 18. Snv and indel calling report generation ( rtg vcfstats )

if calls_report:

    if "calls_report.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of calls_report.log in logs is telling you that the calls report was already produced, please remove calls_report.log if you wish to perform this stage anyway\n")

    
        os.system(
            "%srtg vcfstats %s > %sreport/%s_vcfstats.log" %
            (path_rtg, out, variant_results_file, out, sample_name))

                
                
       
       
        os.system("touch  %slogs/calls_report.log" % (out))


# 19. Html report generation ( Multiqc )


if alignment_report or calls_report or sequencing_report:

    if "multiqc.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of multiqc.log in logs is telling you that the multiqc report was already produced, please remove multiqc.log if you wish to perform this stage anyway\n")

    else:

        os.system("%smultiqc %sreports" % (path_multiqc, out))

        os.system("touch  %slogs/multiqc.log" % (out))

# 20. Annotated variants report generation


if results_report:
	
  if "annovar.log" not in os.listdir(out + "logs"):

        print("WARNING: The absence of annovar.log in logs is telling you that annotation was not peformed, please perform annotation using the -annotation flag if you wish to perform this stage \n")
           
  else:

    if "results_report.log" in os.listdir(out + "logs"):

        print("WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n")

    else:
	
        os.system(
            "zcat %s > %stemp.vcf" %
        
            (variant_results_file, out))

        file = open('%stemp.vcf' % (out), 'r')

        file_lines = file.readlines()

        a = {}

        b = {}

        c = {}

        os.system(
            "cat %stemp.vcf | grep  nonsynonymous > %s/nonsynonymous.vcf" %
            (out, out))

        os.system(
            "cat %stemp.vcf | grep CLINSIG=pathogenic > %s/deleterious.vcf" %
            (out, out))

        file_nonsyn = open('%snonsynonymous.vcf' % (out), 'r')

        file_nonsyn_lines = file_nonsyn.readlines()

        file_del = open('%sdeleterious.vcf' % (out), 'r')

        file_del_lines = file_del.readlines()

        for i in file_lines:

            check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

            if check:

                a[i.split('Gene.refGene=')[1].split(';')[0]] = []

                b[i.split('Gene.refGene=')[1].split(';')[0]] = []

                c[i.split('Gene.refGene=')[1].split(';')[0]] = []

        for i in file_lines:

            check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

            if check:

                a[i.split('Gene.refGene=')[1].split(';')[0]].append(i)

        for i in file_nonsyn_lines:

            check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

            if check:

                b[i.split('Gene.refGene=')[1].split(';')[0]].append(i)

        for i in file_del_lines:

            check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

            if check:

                c[i.split('Gene.refGene=')[1].split(';')[0]].append(i)

        out_file = open('%sreports/output.txt' % (out), 'w')

        for i in a.keys():

            out_file.write(
                'Gene: %s has %s variants, %s nonsynonymous and %s reported to be deleterious\n' %
                (i, len(
                    a[i]), len(
                    b[i]), len(
                    c[i])))

        out_file.close()

        out_file_all = open('%sreports/output_all_variants.txt' % (out), 'w')

        vcf = open('%stemp.vcf' % (out), 'r')

        vcf_lines = vcf.readlines()

        for i in a.keys():

            for j in vcf_lines:

                check1 = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', j, flags=0)

                check = re.search("=%s;" % (i), j, flags=0)

                if check and check1:

                    out_file_all.write(
                        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                        (i,
                         j.split('\t')[0],
                            j.split('\t')[1],
                            j.split('\t')[3],
                            j.split('\t')[4],
                            j.split('Gene.refGene=')[1].split(';')[0],
                            j.split('ExonicFunc.refGene=')[1].split(';')[0],
                            j.split('AAChange.refGene=')[1].split(';')[0],
                            j.split('Func.refGene=')[1].split(';')[0],
                            j.split('CADD_phred=')[1].split(';')[0],
                            j.split('CLINSIG=')[1].split(';')[0],
                            j.split('CLNDSDBID=')[1].split(';')[0]))

        out_file_all.close()

        os.system("touch  %slogs/results_report.log" % (out))




# 21. Starting iobio services

if iobio:
  if "iobio.log" in os.listdir(out + "logs"):

      print("WARNING: The presence of iobio.log in logs is telling you that the iobio services were already started, please remove iobio.log if you wish to start them again\n")

  else:	
  

    os.system(
        "cd %s ; nohup python3 -m http.server %s  >/dev/null 2>&1 &" %
        (path_iobio, port_num))

    os.system("touch  %slogs/iobio.log" % (out))
	
    print("\n\nIobio serces have been started at http://localhost:%s\n\nCopy and paste http://localhost:%s to select the service (vcf, bam, gene) and upload your data into the selected service\n\nIf you want to explore your variant calling results please copy and paste the following URL into your browser and upload the vcf file:\n\n" %(port_num, port_num), end='', flush=True)
    
    if "annovar.log" in os.listdir(out + "logs"):
    	
        print("http://localhost:%s/gene.iobio/?species=Human&rel0=proband&rel1=mother&rel2=father&genes=" %(port_num) , end='', flush=True)
	
        a = {}
	
        os.system(
            "zcat %s > %stemp.vcf" %
        
            (variant_results_file, out))

        file = open('%stemp.vcf' % (out), 'r')

        file_lines = file.readlines()
	
        for i in file_lines:

            check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

            if check:

                a[i.split('Gene.refGene=')[1].split(';')[0]] = []
	
        for i in a.keys() :
    	
            print('%s,' %(i.split(',')[0]), end='', flush=True)
        print(' \n\n')     
    
    else:
		
        print("http://localhost:%s/gene.iobio/" %(port_num))
        
