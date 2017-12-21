#!/bin/bash

reference_file = $1

input1 = $2

input1 = $3

num_cpu = $4

path_to_bwa = $5

path_to_samtools = $6

path_to_java = $7

path_to_picard_jar = $8

path_to_freebayes = $9

path_to_tabix = ${10}

path_to_bcftools = ${11} 


$path_to_bwa/bwa mem -M  -t $num_cpu $reference_file  $input1 $input2 | $path_to_samtools/samtools view  -@ $num_cpu -Shb - | $path_to_samtools/samtools sort -@ $num_cpu -o sorted_reads.bam -

$path_to_java/java -jar  $path_to_picard_jar/picard.jar MarkDuplicates \ 
    INPUT=sorted_reads.bam \ 
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt
	
$path_to_samtools/samtools index dedup_reads.bam

$path_to_freebayes/freebayes -f reference.f dedup_reads.bam > var.vcf

bgzip var.vcf

$path_to_tabix/tabix index -p vcf var.vcf.gz

$path_to_bcftools/bcftools filter -i '' var.vcf.gz | bgzip -c > test_sorted_filtered.vcf.gz ; tabix -fp vcf test_sorted_filtered.vcf.gz