import argparse

import sys, string, os, os.path, re

import paths




path_hisat=paths.path_hisat

path_hisat_index=paths.path_hisat_index

path_samtools=paths.path_samtools

path_freebayes=paths.path_freebayes

path_annovar=paths.path_annovar

path_annovar_db=paths.path_annovar_db

num_cpu=paths.num_cpu

path_bed=paths.path_bed

path_expansionHunter=paths.path_expansionHunter

path_reference=paths.path_reference

path_expansionHunter_jsons=paths.path_expansionHunter_jsons

path_manta=paths.path_manta

path_bedtools=paths.path_bedtools

path_tabix=paths.path_tabix

path_virus_index=paths.path_virus_index

path_bacteria_index=paths.path_virus_index




parser = argparse.ArgumentParser(description='help message here')

parser.add_argument('-format', action="store", dest="format", help='options are bam, sam, fastq, vcf')

parser.add_argument('-paired', action="store", dest="paired", default=1, help='options are 1 for paired end reads and 0 for single end reads')

parser.add_argument('-reference', action="store", dest="reference", help='options are hg19, hg38')

parser.add_argument('-vcf', action="store", dest="vcf", help='complementary vcf file')

parser.add_argument('-in', action="store", dest="input_file", help='input file')

parser.add_argument('-in2', action="store", dest="input_file2", help='for paired end reads only')


parser.add_argument('-out', action="store", dest="out", help='output folder')

parser.add_argument('-expansion', action="store_true", dest="expansion", help='look for expansions? True or False? ', default=False) 

parser.add_argument('-SV', action="store_true", dest="SV", help='look for SVs? True or False? ', default=False) 

parser.add_argument('-BED', action="store_true", dest="BED", help='restrict the analysis to the regione in the bed file? True or False? ', default=False) 

parser.add_argument('-virus', action="store_true", dest="virus", help='perform viral scanning? True or False? ', default=False) 

parser.add_argument('-bacteria', action="store_true", dest="bacteria", help='perform bacteria scanning? True or False? ', default=False) 








args = parser.parse_args()

format=args.format

paired=args.paired

reference=args.reference

input_file=args.input_file

input_file2=args.input_file2

expansion=args.expansion

out=args.out

vcf=args.vcf

SV=args.SV

BED=args.BED

bacteria=args.bacteria

virus=args.virus

# freebayes does support multithreading so I am planning to make it parallel throw the use of python libraries (not yet done). To this purpose I split the bed file into n bed files
if BED==True:
    
    os.system("%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %(path_bedtools, num_cpu,path_bed,out))
    
    
        
if format:
    if format=="fastq" and "freebayes.vcf" not in os.listdir(out):
        
        if paired == 1 :
            
            os.system("%shisat -p %s -x %s -1 %s -2 %s | %ssamtools view -Sb -  |%ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %(path_hisat,num_cpu,path_hisat_index,input_file,input_file2,path_samtools,path_samtools,out,out,path_samtools,out))

            input_file = "%ssorted.bam" %(out)
            
        if paired == "0" :
            
            os.system("%shisat -p %s -x %s -U %s | %ssamtools view -Sb -  |%ssamtools sort -T %stemp.sorted -o %ssorted.bam ; %ssamtools index %ssorted.bam" %(path_hisat,num_cpu,path_hisat_index,input_file,path_samtools,path_samtools,out,out,path_samtools,out))

            input_file = "%ssorted.bam" %(out)
            
    if format=="sam" and "freebayes.vcf" not in os.listdir(out):
        
        
            os.system("%ssamtools view -Sb %s  > %ssorted.bam" %(path_samtools,input_file,out))
        
            input_file = "%ssorted.bam" %(out)




if "freebayes.vcf" not in os.listdir(out):
    
    if vcf :
        
        print ("\nUsing input vcf as variant file\n")

        
    else:
        
        if BED == True:
            
            os.system("%sfreebayes -t %s -f %s %s > %s/freebayes.vcf" %(path_freebayes,path_bed,path_reference,input_file,out))
            
        else:
            
            os.system("%sfreebayes -f %s %s > %s/freebayes.vcf" %(path_freebayes,path_reference,input_file,out))
        
        vcf = "freebayes.vcf"


if "annovar.vcf.hg19_multianno.vcf" not in os.listdir(out):   
    os.system("%stable_annovar.pl --vcfinput %s/%s %s -buildver %s -remove -protocol refGene,dbnsfp30a,clinvar_20170130,avsnp147,cadd -operation g,f,f,f,f -nastring . --outfile %s/annovar.vcf" %(path_annovar,out,vcf,path_annovar_db,reference,out))
    
    
    
if expansion == True and "EH_vcf" not in os.listdir(out):
    
    os.system("%sExpansionHunter --bam %s --ref-fasta %s --repeat-specs %s --vcf %s/EH_vcf --json %s/test.json --log %s/test_log" %(path_expansionHunter,input_file,path_reference,path_expansionHunter_jsons,out,out,out))
    
if SV == True:
    
    if BED == True:
    
        os.system("bgzip -c %s  > %s/temp.bed.gz " %(path_bed,out))
    
        os.system("%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz " %(path_bedtools,out,out))
    
        os.system("%stabix -p bed %s/sorted.bed.gz" %(path_tabix,out))
    
        os.system("%sconfigManta.py --bam %s --referenceFasta %s --runDir %s --callRegions %s/sorted.bed.gz" %(path_manta,input_file,path_reference,out,out))
        
    else:
        
        os.system("%sconfigManta.py --bam %s --referenceFasta %s --runDir %s " %(path_manta,input_file,path_reference,out))

    
    os.system("%srunWorkflow.py -j %s -m local" %(out,num_cpu))
    
    
    
    
    
if virus == True or bacteria == True:
    
    os.system("%ssamtools view -f 4 %s | %ssamtools bam2fq - | gzip - > %suniligned_reads.fastq.gz" %(path_samtools,input_file,path_samtools,out))
    
    if virus == True :
        
        os.system("%shisat -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_virus.bam -" %(path_hisat,num_cpu,path_virus_index,out,path_samtools,path_samtools,out,out))
        
        os.system("%ssamtools index %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt" %(path_samtools,out,path_samtools,out))
        
        os.system("awk \'{printf $1}\' %svirus_stats.txt > %svirus_list.txt ; for i in $(cat %svirus_list.txt); do printf \"$i \"; %ssamtools depth %soutput_virus.bam| head | grep $i | awk \'$3>0 {print $0}\'| wc -l; done > %svirus_coverage_stats.txt" %(out,out,out,path_samtools,out,out))
        
    
    if bateria == True:
    
        os.system("%shisat -p %s -x %s -U %suniligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -" %(path_hisat,num_cpu,path_bacteria_index,out,path_samtools,path_samtools,out,out))

        os.system("%ssamtools index %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt" %(path_samtools,out,path_samtools,out))

        os.system("awk \'{printf $1}\' %sbacteria_stats.txt > %sbacteria_list.txt ; for i in $(cat %sbacteria_list.txt); do printf \"$i \"; %ssamtools depth %soutput_bacteria.bam| head | grep $i | awk \'$3>0 {print $0}\'| wc -l; done > %sbacteria_coverage_stats.txt" %(out,out,out,path_samtools,out,out))





#the following part generates the patient readable output

file=open('%s/annovar.vcf.hg19_multianno.vcf'  %(out), 'r')

file_lines=file.readlines()

a={}
b={}
c={}

os.system("cat %s/annovar.vcf.hg19_multianno.vcf | grep  nonsynonymous > %s/nonsynonymous.vcf" %(out,out))
os.system("cat %s/annovar.vcf.hg19_multianno.vcf | grep CLINSIG=pathogenic > %s/deleterious.vcf" %(out,out))

file_nonsyn=open('%s/nonsynonymous.vcf' %(out),'r')

file_nonsyn_lines=file_nonsyn.readlines()

file_del=open('%s//deleterious.vcf' %(out),'r')

file_del_lines=file_del.readlines()


for i in file_lines:
     check=re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)',i,flags=0)
     if check:
             a[i.split(';')[44].split('=')[1]]=[]
             b[i.split(';')[44].split('=')[1]]=[]
             c[i.split(';')[44].split('=')[1]]=[]
             
             
             
for i in file_lines:
     check=re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)',i,flags=0)
     if check:
             a[i.split(';')[44].split('=')[1]].append(i)
             
for i in file_nonsyn_lines:
     check=re.search(r'(^chr)|(^[0-9,X,Y,M+\t)',i,flags=0)
     if check:
             b[i.split(';')[44].split('=')[1]].append(i)
            
for i in file_del_lines:
    check=re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)',i,flags=0)
    if check:
             c[i.split(';')[44].split('=')[1]].append(i)
             
             
out_file=open('%s/output.txt' %(out),'w' )
             
for i in a.keys():
    
    out_file.write('Gene: %s has %s variants, %s nonsynonymous and %s reported to be deleterious\n' %(i, len(a[i]),len(b[i]),len(c[i])))
    
out_file.close()
    
    
out_file_all=open('%s/output_all_variants.txt' %(out),'w' )

vcf=open('%s/annovar.vcf.hg19_multianno.vcf' %(out) ,'r' )

vcf_lines=vcf.readlines()




for i in a.keys():
    
    #out_file_all.write('%s:\n' %(i))
    
    for j in vcf_lines:
        check1=re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)',j,flags=0)
        check=re.search("=%s;" %(i),j,flags=0)
        
        if check and check1:
            #print(j+"\n")
            #print j.split(';')[43].split('\\')[0]+"\n"
            out_file_all.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(i,j.split('\t')[0],j.split('\t')[1],j.split('\t')[3],j.split('\t')[4],j.split(';')[82],j.split(';')[83],j.split(';')[84],j.split(';')[40],j.split(';')[43].split('\\')[0],j.split(';')[46],j.split(';')[89]))
        
out_file_all.close()
        
        
    
    
    
    
    
    
    


