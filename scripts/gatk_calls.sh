#!/bin/bash

reference_file = $1

input1 = $2

input1 = $3

num_cpu = $4

path_to_bwa = $5

path_to_java = $6

path_to_picard_jar = $7

path_to_gatk_jar = $8

path_to_tabix = $9

path_to_bcftools = ${10} 

wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms//human_9606_b149_GRCh37p13/VCF/GATK/All_20161121.vcf.gz

wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms//human_9606_b149_GRCh37p13/VCF/GATK/All_20161121.vcf.gz.tbi

${path_to_bwa}bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'  -t $num_cpu $reference_file  $input1 $input2 > aligned_reads.sam

${path_to_java}java -jar ${path_to_picard_jar}picard.jar SortSam \ 
    INPUT=aligned_reads.sam \ 
    OUTPUT=sorted_reads.bam \ 
    SORT_ORDER=coordinate 

${path_to_java}java -jar ${path_to_picard_jar}picard.jar MarkDuplicates \ 
	    INPUT=sorted_reads.bam \ 
	    OUTPUT=dedup_reads.bam \
	    METRICS_FILE=metrics.txt
		
${path_to_java}java -jar ${path_to_picard_jar}picard.jar BuildBamIndex \ 
		    INPUT=dedup_reads.bam
			
${path_to_java}java -jar ${path_to_gatk_jar}GenomeAnalysisTK.jar \ 
			-T BaseRecalibrator \ 
			-R $reference_file \ 
			-I dedup_reads.bam \ 
			-L 20 \ 
			-knownSites All_20161121.vcf.gz 
			-o recal_data.table
		
${path_to_java}java -jar ${path_to_gatk_jar}GenomeAnalysisTK.jar \ 
			-T BaseRecalibrator \ 
			-R $reference_file \ 
			-I dedup_reads.bam \ 
			-L 20 \ 
			-knownSites All_20161121.vcf.gz
			-BQSR recal_data.table \ 
			-o post_recal_data.table 
	
${path_to_java}java -jar ${path_to_gatk_jar}GenomeAnalysisTK.jar \ 
			-T PrintReads \ 
			-R $reference_file \ 
			-I dedup_reads.bam \ 
			-L 20 \ 
			-BQSR post_recal_data.table \ 
			-o recal_reads.bam
			
${path_to_java}java -jar ${path_to_picard_jar}picard.jar BuildBamIndex \ 
					    INPUT=dedup_reads.bam
		
${path_to_java}java -jar ${path_to_gatk_jar}GenomeAnalysisTK.jar \ 
			-T HaplotypeCaller \ 
			-R $reference_file \ 
			-I recal_reads.bam \  
			-L 20 \ 
			--genotyping_mode DISCOVERY \ 
			-o raw_variants.vcf
			
bgzip raw_variants.vcf

${path_to_tabix}tabix -p vcf raw_variants.vcf.gz

${path_to_bcftools}bcftools filter -i 'QUAL > 30 & QUAL/AO > 10' raw_variants.vcf.gz | bgzip -c > test_gatk_sorted_filtered.vcf.gz ; ${path_to_tabix}tabix -p vcf test_gatk_sorted_filtered.vcf.gz
