#!/bin/bash

reference_file=$1

input1=$2

input1=$3

num_cpu=$4

path_to_bwa=$5

path_to_samtools=$6

path_to_java=$7

path_to_picard_jar=$8

path_to_freebayes=$9

path_to_tabix=${10}

path_to_bcftools=${11}


${path_to_bwa}bwa mem -M  -t $num_cpu $reference_file  $input1 $input2 | ${path_to_samtools}samtools view  -@ $num_cpu -Shb - | ${path_to_samtools}samtools sort -@ $num_cpu -o sorted_reads.bam -

${path_to_java}java  -XX:ParallelGCThreads=${num_cpu} -jar  ${path_to_picard_jar}picard.jar MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt

${path_to_samtools}samtools index -@ $num_cpu dedup_reads.bam

#${path_to_freebayes}freebayes-parallel <(${path_to_freebayes}fasta_generate_regions.py ${reference_file}.fai 100000) $num_cpu -f $reference_file dedup_reads.bam > var.vcf

${path_to_freebayes}freebayes -f $reference_file dedup_reads.bam > var.vcf

bgzip var.vcf

${path_to_tabix}tabix -p vcf var.vcf.gz

${path_to_bcftools}bcftools filter -i 'QUAL > 30 & QUAL/AO > 10' var.vcf.gz | bgzip -c > test_sorted_filtered.vcf.gz ; ${path_to_tabix}tabix -p vcf test_sorted_filtered.vcf.gz

