import argparse

import sys, string, os, os.path, re

import paths, subprocess


#Define paths variables from paths.py

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

path_bwa = paths.path_bwa

path_bwa_index = paths.path_bwa_index

path_rtg = paths.path_rtg

path_bcftools = paths.path_bcftools



#Parse options variables from command line

parser = argparse.ArgumentParser( description = 'help message here' )

parser.add_argument( '-filter_string' , action = "store" , dest = "filter_string" , default = "" , help = 'filter string, eg GQ>20 & DP>10' )

parser.add_argument( '-mode' , action = "store" , dest = "mode" , default = "fast" , help = 'options are fast, normal, intensive' )

parser.add_argument( '-format' , action = "store" , dest = "format" , help = 'options are bam, sam, fastq, vcf' )

parser.add_argument( '-paired', action = "store" , dest = "paired" , default = 1 , help = 'options are 1 for paired end reads and 0 for single end reads' )

parser.add_argument( '-reference' , action = "store" , dest = "reference" , help = 'options are hg19, hg38' )

parser.add_argument( '-vcf' , action = "store" , dest = "vcf" , help = 'complementary vcf file' )

parser.add_argument( '-in' , action = "store" , dest = "input_file" , help = 'input file' )

parser.add_argument( '-in2' , action = "store" , dest = "input_file2" , help = 'for paired end reads only' )

parser.add_argument( '-iobio' , action = "store_true" , dest = "iobio" , help = 'wanna start iobio services at the end of the analysis?' , default = False ) 

parser.add_argument( '-alignment', action = "store_true" , dest = "alignment" , help = 'Align your data to a reference genome? True or False? ' , default = False ) 

parser.add_argument( '-out', action = "store" , dest = "out" , help = 'output folder' )

parser.add_argument( '-expansion' , action = "store_true" , dest = "expansion" , help = 'look for expansions? True or False? ' , default = False ) 

parser.add_argument( '-SV' , action = "store_true" , dest = "SV" , help = 'look for SVs? True or False? ' , default = False ) 

parser.add_argument( '-BED' , action = "store_true" , dest = "BED" , help = 'restrict the analysis to the regione in the bed file? True or False? ' , default = False ) 

parser.add_argument( '-virus' , action = "store_true" , dest = "virus" , help = 'perform viral scanning? True or False? ' , default = False ) 

parser.add_argument( '-bacteria' , action = "store_true" , dest = "bacteria" , help = 'perform bacteria scanning? True or False? ' , default = False ) 

parser.add_argument( '-variantcalling', action = "store_true" , dest = "variantcalling" , help = 'wanna perform variant calling? True or False? ' , default = False ) 

parser.add_argument( '-annotation', action = "store_true" , dest = "annotation" , help = 'wanna perform annotation? True or False? ' , default = False ) 

parser.add_argument( '-results_report', action = "store_true" , dest = "results_report" , help = 'wanna perform report? True or False? ' , default = False ) 

parser.add_argument( '-alignment_report', action = "store_true" , dest = "alignment_report" , help = 'wanna perform alignment report? True or False? ' , default = False ) 

parser.add_argument( '-sequencing_report', action = "store_true" , dest = "sequencing_report" , help = 'wanna perform sequencing report? True or False? ' , default = False ) 

parser.add_argument( '-calls_report', action = "store_true" , dest = "calls_report" , help = 'wanna perform calls report? True or False? ' , default = False ) 

parser.add_argument( '-sample_name' , action = "store" , dest = "sample_name" , default = "sample" , help = 'specify sample name' )

parser.add_argument( '-rm_dup', action = "store_true" , dest = "rm_dup" , help = 'wanna remove duplicate? True or False? ' , default = False ) 

parser.add_argument( '-debug' , action = "store_true" , dest = "debug" , help = 'wanna run the pipeline in debug mode?' , default = False ) 


#Define parsed command line options 

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

os.system( "mkdir  %slogs ; mkdir  %sreports ; mkdir  %sresults" %( out , out , out ) )


#If bed is not provided, it generates one starting from the reference genome index. The pipeline uses a bed file to split the analyses in several subprocesses  

if BED == True :
    
    os.system( "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %( path_bedtools , num_cpu , path_bed , out ) )
    
else :
    
    os.system( "cat %s.fai | awk '{print $1\"\t0\t\"$2} > %sreference.bed'" %( path_reference , out) )
    
    path_bed = "%sreference.bed" %(out)
    
    os.system( "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %( path_bedtools , num_cpu , path_bed , out ) )

    BED = True
    
#perform alignment if input sequencing data is in fastq format
        
if alignment == True :        
    
    if format == "fastq" and "alignment.log" not in os.listdir( out + "logs" ) : 
        
        if paired == 1 :
            
            if mode == "fast" :
                
                os.system( "%shisat2 --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %ssamtools view -Sb -  | %ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %( path_hisat , num_cpu , path_hisat_index , input_file , input_file2 , path_samtools , path_samtools , out , out , path_samtools , out ) )

                bam_file = "%ssorted.bam" %( out )
                
                os.system("touch  %salignment.log" %( out ) )
                
            if mode == "normal" or mode == "intensive" :

                os.system( "%shisat2 --rg-id 4 --rg LB:lib1 --rg PL:illumina  --rg PU:unit1 --rg SM:20 --no-softclip --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted.bam ; %ssamtools index %ssorted.bam" %( path_hisat , num_cpu , path_hisat_index , input_file , input_file2 , path_samtools , path_samtools , out , out , path_samtools , out ) )

                os.system( "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq" %( path_samtools , num_cpu , out , out ) )

                os.system( "%sbwa mem -R '@RG\tID:4\tLB:lib1\tPL:illumina\tRGPU:unit1\tSM:20' -t %s %s %sunaligned_reads.fq | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted_bwa.bam ; %ssamtools index %ssorted_bwa.bam " %( path_bwa , num_cpu , path_bwa_index , out , path_samtools , path_samtools , out , out , path_samtools , out ) )

                os.system( "%ssamtools view -H %ssorted.bam > %sheader.txt" %( path_samtools , out , out ) )

                os.system( "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam" %( path_samtools , num_cpu , out , out , out , out ) ) 

                if debug == False :
                    
                        os.system( "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq %sheader.txt" %( out , out , out , out ) )

                os.system( "%ssamtools index %ssorted_merged.bam" %( path_samtools , out ) )

                bam_file = "%ssorted_merged.bam" %( out )

                os.system( "touch  %slogs/alignment.log" %( out ) )
                
    
        
            
                   
        if paired == "0" :
            
            if mode == "fast" :
            
                os.system( "%shisat2 --no-spliced-alignment--remove-chrname -p %s -x %s -U %s | %ssamtools view -Sb -  | %ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %( path_hisat , num_cpu , path_hisat_index , input_file , path_samtools , path_samtools , out , out , path_samtools , out ) )

                bam_file = "%ssorted.bam" %( out )
                
                os.system( "touch  %salignment.log" %( out ) )
                
                
            if mode == "normal" or mode == "intensive" :  
            
                os.system( "%shisat2 --rg-id 4 --rg LB:lib1 --rg PL:illumina  --rg PU:unit1 --rg SM:20 --no-softclip --no-spliced-alignment -p %s -x %s -U %s | %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted.bam ; %ssamtools index %ssorted.bam" %( path_hisat , num_cpu , path_hisat_index , input_file , path_samtools , path_samtools , out , out , path_samtools , out ) )
 
                os.system( "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq" %( path_samtools , num_cpu , out , out ) )

                os.system( "%sbwa mem -R '@RG\tID:4\tLB:lib1\tPL:illumina\tRGPU:unit1\tSM:20' -t %s %s %sunaligned_reads.fq| %ssamtools view -Sb -  | %ssamtools sort -T %s -o %ssorted_bwa.bam ; %ssamtools index %ssorted_bwa.bam " %( path_bwa , num_cpu , path_bwa_index , out , path_samtools , path_samtools , out , out , path_samtools , out ) )

                os.system( "%ssamtools view -H %ssorted.bam > %sheader.txt" %( path_samtools , out , out) )

                os.system( "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam" %( path_samtools , num_cpu , out , out , out , out ) ) 

                os.system( "%ssamtools index %ssorted_merged.bam " %( path_samtools , out ) )
                
                if debug == False :
                    
                        os.system( "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq " %( out , out , out ) )

                bam_file = "%ssorted_merged.bam" %( out )

                os.system( "touch  %slogs/alignment.log" %( out ) )
                
    else:
        
        if format != "fastq" :
            
            print( "WARNING: Fastq format input data is requiered if you want to perform the alignment stage\n" )
            
        if "alignment.log" in os.listdir( out + "logs" ) :
            
            print( "WARNING: The presence of alignment.log in logs is telling you that the alignment was already peformed, please remove alignment.log if you wish to perform this stage anyway\n" )
    
            if mode == "normal" or mode == "intensive" : 
            
                bam_file = "%ssorted_merged.bam" %( out )
            
            if mode == "fast" :
                
                bam_file = "%ssorted.bam" %( out )
                
#Remove duplicates       
                
if rm_dup == True :
    
    
    if "rm_dup.log" in os.listdir( out + "logs" ) :
        
        print( "WARNING: The presence of rm_dup.log in logs is telling you that the duplicate removal was already peformed, please remove rm_dup.log if you wish to perform this stage anyway\n" )
    
    
    else : 
    
        if paired == 0 :
        
            single_end_option = "-S"
        
        else :
        
            single_end = ""
        
        os.system( "%ssamtools rmdup %s %s %ssorted_merged_rmdup.bam" %(path_samtools , single_end_option , bam_file , out ) )
    
        bam_file = "%ssorted_merged_rmdup.bam" %( out )
    
        if debug == False :
        
            os.system( "rm %ssorted_merged.bam" %( out ) )
            
        s.system( "touch  %slogs/rm_dup.log" %( out ) ) 
                
#Convert input sam file into bam 
            
if format == "sam" and "sam2bam.log" not in os.listdir( out + "logs" ) :
        
    os.system( "%ssamtools view -Sb %s  > %ssorted.bam" %( path_samtools , input_file , out ) )
        
    bam_file = "%ssorted.bam" %( out )
    
    os.system( "touch  %slogs/sam2bam.log" %( out ) ) 


#Perform variant calling (snps and indels)

if variantcalling == True:
    
    if vcf :
        
        print( "WARNING: Using input vcf as variant file. Do not provide vcf file if you wish to perform variant calling\n" )
        
    else :
        
        if "VC.log" in os.listdir( out + "logs" ) :
            
            print( "WARNING: The presence of VC.log in logs is telling you that the variant calling was already peformed, please remove VC.log if you wish to perform this stage anyway\n" )
            
            vcf = "%s%s_sorted.vcf.gz" %( out , sample_name )
            
        else:
            
            if mode == "intensive"  and "VC_gatk.log" not in os.listdir( out + "logs" ) : 
    
                counter = 1
    
                ps = []
    
                while counter < int( num_cpu ) +1 :   
                    
                    command = "%ssamtools mpileup --max-depth 10000 -AQ 0 -l %stemp%s.bed %s | awk '$5 ~/[\+,\-]/ {print $1\"\t\"$2-1\"\t\"$2}' > %smpileup_positions%s.bed" %( path_samtools , out , str( counter ) , bam_file ,  out , str( counter ) )
    
                    proc_mpileup = subprocess.Popen( command , shell = True )
    
                    ps.append( proc_mpileup )           
    
                    counter += 1
                
                for proc_mpileup in ps :
                
                    proc_mpileup.wait()
                    
                counter = 1
    
                ps = []
    
                while counter < int( num_cpu ) +1 :
    
                    command = "%sjava -jar %sGenomeAnalysisTK.jar -R %s -T HaplotypeCaller -I %s -L %smpileup_positions%s.bed -o %sgatk_indels%s.vcf" %( path_java , path_gatk , path_reference , bam_file , out , str( counter ) , out , str( counter ) )
    
                    proc_gatk = subprocess.Popen( command , shell = True )
    
                    ps.append( proc_gatk )           
    
                    counter += 1
                
                for proc_gatk in ps :
                
                    proc_gatk.wait()
                
                os.system( "cat %sgatk_indels1.vcf | grep \"^#\" >> %sgatk_indels_merged.vcf ; for i in $(ls %s | grep gatk_indels | grep -v idx); do cat %s$i | grep -v \"^#\" >> %sgatk_indels_merged.vcf; done" %( out , out , out , out , out ) )
    
                os.system( "%sbedtools sort -header -i %sgatk_indels_merged.vcf > %sgatk_indels_sorted_merged.vcf" %( path_bedtools , out , out ) )
    
                os.system( "%svcftools  --vcf %sgatk_indels_sorted_merged.vcf --minGQ 20 --minDP 2 --keep-only-indels --recode --recode-INFO-all --out %sindels_only" %( path_vcftools , out , out ) )
                    
                if debug == False :
        
                    os.system( "rm %sgatk_indels* %smpileup_positions* %sindels_only.log" %( out , out , out ) )
                
                os.system( "touch  %slogs/VC_gatk.log" %( out ) )
                 
                    
            if "VC_freebayes.log" not in os.listdir( out + "logs" ) : 
    
                counter = 1
                
                ps = []
            
                os.system( "date" )
            
                while counter < int(num_cpu)+1:

                    command = path_freebayes + "freebayes" + " -t " + out + "temp" + str( counter ) + ".bed" + " -f " + path_reference + " -b " + bam_file + " > " + out + "/freebayes" + str( counter ) + ".vcf"

                    proc_freebayes = subprocess.Popen( command , shell = True )
	   
                    ps.append( proc_freebayes )           
               
                    counter += 1
               
               
                for proc_freebayes in ps :
                
                    proc_freebayes.wait()
            
                os.system( "date" )
            
                os.system( "cat %sfreebayes1.vcf | grep \"^#\" >> %smerged.vcf ; for i in $(ls %s | grep freebayes) ; do cat %s$i | grep -v \"^#\" >> %smerged.vcf; done" %( out , out , out , out , out ) )
            
                os.system( "%sbedtools sort -header -i %smerged.vcf > %sfreebayes.vcf" %( path_bedtools , out , out ) )
        
        
                vcf = "freebayes.vcf"
            
            
        
                if mode == "intensive":
        
                    os.system( "%svcftools  --vcf %sfreebayes.vcf --minGQ 20 --minDP 2 --remove-indels --recode --recode-INFO-all --out %sSNPs_only" %( path_vcftools , out , out ) )
                    
                    os.system( "bgzip  %sSNPs_only.recode.vcf ; bgzip %sindels_only.recode.vcf " %(  out , out ) )
            
                    os.system( "%stabix -p vcf %sSNPs_only.recode.vcf.gz ; %stabix -p vcf %sindels_only.recode.vcf.gz " %( path_tabix , out , path_tabix , out ) )
                    
                           
                else:
                    
                      os.system( "cat %sfreebayes.vcf | bgzip -c > %s%s_sorted.vcf.gz" %( out , out , sample_name ) )
                      
                      os.system( "tabix -p vcf %s%s_sorted.vcf.gz" %( out , sample_name ) )
            
                
                
                if debug == False :
        
                    os.system( "rm %sfreebayes* %sSNPs_only.log %stemp* %smerged.vcf" %( out , out , out , out) )
                
                os.system( "touch  %slogs/VC_freebayes.log" %( out ) )
            
           
                
            os.system( "touch  %slogs/VC.log" %( out ) )
            
            
if filter_string :
    
    if mode == "intensive" :
    
        os.system("%sbcftools filter -i \'%s\' %sSNPs_only.recode.vcf.gz | bgzip -c > %stemp1.vcf.gz ; mv %stemp1.vcf.gz %sSNPs_only.recode.vcf.gz ; %tabix -fp vcf %sSNPs_only.recode.vcf.gz" (path_bcftools , filter_string , out , out , out , out , path_tabix , out ) )
            
        os.system("%sbcftools filter -i \'%s\' %sindels_only.recode.vcf.gz | bgzip -c > %stemp2.vcf.gz ; mv %stemp2.vcf.gz %sindels_only.recode.vcf.gz ; %tabix -fp vcf %sindels_only.recode.vcf.gz" (path_bcftools , filter_string , out , out , out , out ,  path_tabix , out ) )
    
    else :
       
       os.system("%sbcftools filter -i \'%s\' %s%s_sorted.vcf.gz| bgzip -c > %stemp1.vcf.gz ; mv %stemp1.vcf.gz %s%s_sorted.vcf.gz ; %tabix -fp vcf %s%s_sorted.vcf.gz" (path_bcftools , filter_string , out , sample_name , out , out , sample_name, out , path_tabix , out , sample_name ) )
        
                
if annotation == True :
    
    if "annovar.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of annovar.log in logs is telling you that annotation was already peformed, please remove annovar.log if you wish to perform this stage anyway\n" )
        
        
    else:  
        
        if mode == "intensive":
    
            os.system( "%sjava -jar %sGenomeAnalysisTK.jar -T CombineVariants -R %s --variant %sSNPs_only.recode.vcf.gz --variant %sindels_only.recode.vcf.gz -o  %s%s.vcf --genotypemergeoption UNSORTED" %( path_java , path_gatk , path_reference, out , out , out , sample_name ) )
        
            os.system( "%svcf-sort  %s%s.vcf | bgzip -c > %s%s_sorted.vcf.gz" %( path_vcftools , out , sample_name , out , sample_name ) )
            
            if debug == False :
    
                os.system( "rm %ssample.vcf* " %( out ) )
                    
            os.system( "%stabix -p vcf %s%s_sorted.vcf.gz" %( path_tabix , out , sample_name ) )
        
        os.system( "%stable_annovar.pl --vcfinput %s%s_sorted.vcf.gz %s -buildver %s -remove -protocol refGene,dbnsfp30a,clinvar_20170130,avsnp147,cadd -operation g,f,f,f,f -nastring . --outfile %s/annovar.vcf" %( path_annovar , out , sample_name , path_annovar_db , reference , out ) )
            
        os.system("rm %sannovar.vcf.hg19_multianno.txt %sannovar.vcf.avinput" %( out , out ))
        
        os.system( "touch  %slogs/annovar.log" %( out ) )
    
if expansion == True :
    
    
    if "EH.log" in os.listdir( out ) :
        
        print( "WARNING: The presence of EH.log in logs is telling you that the expansion scan was already peformed, please remove SV.log if you wish to perform this stage anyway\n" )
        
    else :
        
        os.system( "%sExpansionHunter --bam %s --ref-fasta %s --repeat-specs %s --vcf %s/EH_vcf --json %s/test.json --log %s/test_log" %( path_expansionHunter , bam_file , path_reference , path_expansionHunter_jsons , out , out , out ) )
    
        os.system( "mv %s/EH_vcf %s/results/%s_expansions.vcf ; bgzip %s/results/%s_expansions.vcf ; %stabix -p vcf %s/results/%s_expansions.vcf.gz" %( out , out , sample_name , out , sample_name , path_tabix , out , sample_name ) )
        
        os.system( "touch  %slogs/EH.log" %( out ) )
    
if SV == True:
    
    if "SV.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of SV.log in logs is telling you that structural variant calling was already peformed, please remove SV.log if you wish to perform this stage anyway\n" )
    
    
    else:
        
        if BED == True:
    
            os.system( "bgzip -c %s  > %s/temp.bed.gz" %( path_bed , out ) )
    
            os.system( "%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz" %( path_bedtools , out , out ) )
            
            os.system( "%stabix -p bed %s/sorted.bed.gz" %( path_tabix , out ) )
            
            os.system( "mkdir %smanta" %( out ) )
    
            os.system( "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta --callRegions %s/sorted.bed.gz" %( path_manta , bam_file , path_reference , out , out ) )
        
        else:
            
            os.system( "mkdir %smanta" %( out ) )
            
            os.system( "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta" %( path_manta , bam_file , path_reference , out ) )

        
        os.system( "%smanta/runWorkflow.py -j %s -m local" %( out , num_cpu ) )
        
        os.system( "mv %s/manta/results/variants/diploidSV.vcf.gz  %s/results/%s_SV.vcf.gz" %( out , out , sample_name ) )
        
        os.system( "mv %s/manta/results/variants/diploidSV.vcf.gz.tbi  %s/results/%s_SV.vcf.gz.tbi" %( out , out , sample_name) )
    
        if debug == False :

            os.system( "rm -r %stemp.bed.gz  %ssorted.bed.gz %smanta" %( out , out , out ) )
            
            
        
        os.system( "touch  %slogs/SV.log" %( out ) )
    
    
    
if virus == True or bacteria == True :
    
    
    if "microbes.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of microbes.log in logs is telling you that microbes scanning was already peformed, please remove microbes.log if you wish to perform this stage anyway\n" )
    
    else:
        
        os.system( "%ssamtools view -hf 4 %s | %ssamtools bam2fq - | gzip -c - > %suniligned_reads.fastq.gz" %( path_samtools , bam_file , path_samtools , out ) )
    
        if virus == True :
        
            os.system( "%shisat2 --no-spliced-alignment -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_virus.bam -" %( path_hisat , num_cpu , path_virus_index , out , path_samtools , path_samtools , out , out ) )
        
            os.system( "%ssamtools index %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt" %( path_samtools , out , path_samtools , out , out ) )
        
            os.system( "awk \'{print $1}\' %svirus_stats.txt > %svirus_list.txt ; for i in $(cat %svirus_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_virus.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %svirus_coverage_stats.txt" %( out , out , out , path_samtools , out , out ) )
        
            virus_coverage_stats = open( '%svirus_coverage_stats.txt' %( out ) , 'r' )

            virus_coverage_stats_lines = virus_coverage_stats.readlines()
            
            virus_stats_file = open( '%svirus_stats.txt' %( out ) , 'r' )

            virus_stats_file_lines = virus_stats_file.readlines()
            
            i = 0
            
            virus_results = open('%sresults/virus_results.txt' %( out ) , 'w' ) 
            
            virus_results.write( "Id\tGenome_lenth\tNumber_of_reds\tCoverage\n" )
            
            while i < len(virus_coverage_stats_lines) : 
                
                virus_results.write( "%s\t%s\t%s\t%s\n" %( virus_stats_file_lines[i].split('|')[1] ,virus_stats_file_lines[i].split('\t')[1] , virus_stats_file_lines[i].split('\t')[2] , virus_coverage_stats_lines[i].split(' ')[1].strip() ) )
            
                i += 1
                
            virus_results.close()
                
            if debug == False :
    
                os.system( "rm %svirus_stats.txt  %svirus_coverage_stats.txt" %( out , out  ) )
            
            
        if bacteria == True :
    
            os.system( "%shisat2 --no-spliced-alignment -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -" %( path_hisat , num_cpu , path_bacteria_index , out , path_samtools , path_samtools , out , out ) )

            os.system( "%ssamtools index %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt" %( path_samtools , out , path_samtools , out , out ) )

            os.system( "awk \'{print $1}\' %sbacteria_stats.txt > %sbacteria_list.txt ; for i in $(cat %sbacteria_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_bacteria.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %sbacteria_coverage_stats.txt" %( out , out , out , path_samtools , out , out ) )

            bacteria_coverage_stats = open( '%sbacteria_coverage_stats.txt' %( out ) , 'r' )

            bacteria_coverage_stats_lines = bacteria_coverage_stats.readlines()
            
            bacteria_stats_file = open( '%sbacteria_stats.txt' %( out ) , 'r' )

            bacteria_stats_file_lines = bacteria_stats_file.readlines()
            
            i = 0
            
            bacteria_results = open('%sresults/bacteria_results.txt' %( out ) , 'w' )
            
            bacteria_results.write( "Id\tGenome_lenth\tNumber_of_reds\tCoverage\n" )
            
            while i < len(bacteria_coverage_stats_lines) : 
                
                bacteria_results.write( "%s\t%s\t%s\t%s\n" %( bacteria_stats_file_lines[i].split('|')[1] , bacteria_stats_file_lines[i].split('\t')[1] , bacteria_stats_file_lines[i].split('\t')[2] , bacteria_coverage_stats_lines[i].split(' ')[1].strip() ) )
            
                i += 1

            bacteria_results.close()
            
            if debug == False :
    
                os.system( "rm %sbacteria_stats.txt  %sbacteria_coverage_stats.txt" %( out , out  ) )
            

        os.system( "touch  %slogs/microbes.log" %( out ) )


if alignment_report == True :
    
    if "alignment_report.log" in os.listdir( out + "logs" ) : 

        print( "WARNING: The presence of alignment_report.log in logs is telling you that the alignment report was already produced, please remove alignment_report.log if you wish to perform this stage anyway\n" )
        
        
    else :
    
        os.system( "%ssamtools flagstat -@ %s %s > %sreport/%s_flagstat.log" %( path_samtools , num_cpu , bam_file , out, sample_name) )
        
        os.system( "%ssamtools stats -@ %s %s > %sreport/%s_stats.log" %( path_samtools , num_cpu , bam_file , out, sample_name ) )
        
        os.system( "touch  %slogs/alignment_report.log" %( out ) )

    
if sequencing_report == True :
    
    
    if "sequencing_report.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of sequencing_report.log in logs is telling you that the sequence report was already produced, please remove sequencing_report.log if you wish to perform this stage anyway\n" )
        
        
    else :
    
        if  path_java != "" :
    
            java_option = "-j " + path_java + " "
    
        os.system( "%sfastqc %s -o %sreports -f %s -t %s %s" %( path_fastqc , java_option , out , format , num_cpu , nput_file ) )
    
        os.system( "touch  %slogs/sequencing_report.log" %( out ) )
        
    
    
if calls_report == True :
    
    if "calls_report.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of calls_report.log in logs is telling you that the calls report was already produced, please remove calls_report.log if you wish to perform this stage anyway\n" )
        
    
    else :    
        
        os.system( "%srtg vcfstats %s%s_sorted.vcf.gz > %sreport/%s_vcfstats.log" %( path_rtg , out , sample_name , out, sample_name ) )
    
        os.system( "touch  %slogs/calls_report.log" %( out ) )
        
        
if alignment_report == True or calls_report == True or sequencing_report == True :
    
    if "multiqc.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of multiqc.log in logs is telling you that the multiqc report was already produced, please remove multiqc.log if you wish to perform this stage anyway\n" )
        
    
    else :
    
        os.system( "%smultiqc %sreports" %( path_multiqc , out ) )
    
        os.system( "touch  %slogs/multiqc.log" %( out ) )
        

if results_report == True :
    
    
    if "results_report.log" in os.listdir( out + "logs" ) : 
        
        print( "WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n" )
        
    else :

        file=open( '%s/annovar.vcf.hg19_multianno.vcf'  %( out ) , 'r' )

        file_lines = file.readlines()

        a = {}

        b = {}

        c = {}
 
        os.system( "cat %s/annovar.vcf.hg19_multianno.vcf | grep  nonsynonymous > %s/nonsynonymous.vcf" %( out , out ) )

        os.system( "cat %s/annovar.vcf.hg19_multianno.vcf | grep CLINSIG=pathogenic > %s/deleterious.vcf" %( out , out ) )

        file_nonsyn = open( '%snonsynonymous.vcf' %( out ) ,'r' )

        file_nonsyn_lines = file_nonsyn.readlines()
    
        file_del = open( '%sdeleterious.vcf' %( out ) , 'r' )

        file_del_lines = file_del.readlines()


        for i in file_lines :
    
             check = re.search( r'(^chr)|(^[0-9,X,Y,M]+\t)' , i , flags=0 )
     
             if check :
     
                 a[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ] = []
     
                 b[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ] = []
     
                 c[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ] = []
             
             
             
        for i in file_lines :
    
     
             check = re.search( r'(^chr)|(^[0-9,X,Y,M]+\t)' , i , flags=0 )
         
             if check:
     
                 a[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ].append( i )
                 
        for i in file_nonsyn_lines :
     
             check = re.search( r'(^chr)|(^[0-9,X,Y,M]+\t)' , i , flags = 0 )
     
             if check :
      
                 b[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ].append( i )
            
        for i in file_del_lines :
    
            check = re.search( r'(^chr)|(^[0-9,X,Y,M]+\t)' , i , flags=0 )
    
            if check :
    
                 c[ i.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ] ].append( i )
             
             
        out_file = open( '%sreports/output.txt' %( out ) , 'w' )
             
    
        for i in a.keys() :
    
            out_file.write( 'Gene: %s has %s variants, %s nonsynonymous and %s reported to be deleterious\n' %( i, len( a[ i ] ) , len( b[ i ] ) , len( c[ i ] ) ) )
    
        out_file.close()
    
        out_file_all = open('%sreports/output_all_variants.txt' %( out ) , 'w' )

        vcf = open( '%s/annovar.vcf.hg19_multianno.vcf' %( out ) , 'r' )

        vcf_lines = vcf.readlines()



        for i in a.keys() :
        
            for j in vcf_lines :
        
                check1 = re.search( r'(^chr)|(^[0-9,X,Y,M]+\t)' , j , flags = 0 )
        
                check = re.search( "=%s;" %( i ) , j , flags = 0 )
        
                if check and check1 :

                    out_file_all.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( i , j.split( '\t' )[ 0 ] , j.split( '\t' )[ 1 ] , j.split( '\t' )[ 3 ] , j.split( '\t' )[ 4 ] , j.split( 'Gene.refGene=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'ExonicFunc.refGene=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'AAChange.refGene=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'Func.refGene=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'CADD_phred=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'CLINSIG=' )[ 1 ].split( ';' )[ 0 ]  , j.split( 'CLNDSDBID=' )[ 1 ].split( ';' )[ 0 ]  ) )
        
        out_file_all.close()
        
        
        os.system( "touch  %slogs/results_report.log" %( out ) )
        
os.system( "mv %s/annovar.vcf.hg19_multianno.vcf %sresults/sample_annotated.vcf ; bgzip %sresults/sample_annotated.vcf ; %stabix -fp vcf %sresults/sample_annotated.vcf.gz" %( out , out , out , path_tabix , out ) )
   
   
        
if iobio == True:
        
        os.system( "pushd %s ; nohup python3.5 -m SimpleHTTPServer %s  &" %( path_iobio , port_num) )
