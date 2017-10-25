################################################################
# Program: DNAscan - script to run DNAscan on a list of samples
# Version 1.0
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

################################################################
# Script structure:
# 1. Import python modules
# 2. Define paths viriables
# 3. Define options from command line
# 4. Parse options from command line
# 5. Run DNAscan for each line in the input sample list
#   5.1 Create create DNAscan input file option string per line in the input list
#   5.2 Create working dir tree
#   5.3 Run DNAscan for one sample
################################################################

import argparse , os  ,  paths , os.path 

from argparse import RawTextHelpFormatter

# 2. Define paths viriables

python_path = paths.python_path

dnascan_dir = paths.dnascan_dir


# 3. Define options from command line

parser = argparse.ArgumentParser(prog='python analyse_list_of_samples.py ', usage='%(prog)s -format "string" -paired "string" -sample_list "string" -out_dir "string" -option_string "string"', description = '############Help Message############ \n\nThis is a script to run DNAscan on a list of samples. Each line of the list must contain the path to one sample. If samples are in paired reads in fastq format and you have two files per sample, these will have to be on the same line spaced bt a tab.\n\n E.g. sample.1.fq.gz  sample.2.fq.gz\n\nDNAscan uses the file paths.py to locate the needed tools and files. Please make sure your paths.py file is properly filled \n\nUsage example: \n\npython alalyse_list_of_files.py -option_string "-format fastq -mode intensive -reference hg19 -alignment -variantcalling -annotation" -out_dir /path/to/dir -sample_list list.txt -format bam\n\nPlease check the following list of required options\n\n################################################', formatter_class=RawTextHelpFormatter)

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument( '-option_string' , required=True , action = "store" , dest = "option_string" , default = "" , help = 'string of option to pass to the main script. Do not include neither -out option nor -in option . E.g. "-format fastq -mode intensive -reference hg19 -alignment -variantcalling -annotation"' )

requiredNamed.add_argument( '-out_dir' , required=True , action = "store" , dest = "out_dir" , default = "" , help = 'working directory [string]' )

requiredNamed.add_argument( '-sample_list' , required=True , action = "store" , dest = "sample_list" , default = "" , help = 'file containing the list of samples to Analyse [string]' )

requiredNamed.add_argument( '-format' , required=True , action = "store" , dest = "format" , help = 'options are bam, sam, fastq, vcf [string]' )

requiredNamed.add_argument( '-paired' , required=True , action = "store" , dest = "paired" , default = "1" , help = 'options are 1 for paired end reads and 0 for single end reads [string]' )


# 4. Parse options from command line

args = parser.parse_args()

option_string = args.option_string

out_dir = args.out_dir

sample_list = args.sample_list

format = args.format

paired = args.paired

list_file = open( "%s" %(sample_list) , 'r' )

list_file_lines = list_file.readlines()



# 5. Run DNAscan for each line in the input sample list

for sample in list_file_lines :
    
    # 5.1 Create create DNAscan input file option string per line in the input list
    
    if paired == "1" and format == "fastq" :
        
        input_file_string = "-in %s -in2 %s" %(sample.split('\t')[0] , sample.split('\t')[1].strip())
                
        sample_name = sample.split('\t')[0].split("/")[-1].split("1.f")[-2]
        
    else :
        
           input_file_string = "-in %s" %( sample.strip() ) 
           
           sample_name = sample.split('.')[-2].split("/")[-1] 
           
    # 5.2 Create working dir tree
    
    os.system( "mkdir %s ; mkdir %s/%s" %( out_dir , out_dir , sample_name ) )
    
    # 5.3 Run DNAscan for one sample
    
    os.system( "%s %sDNAscan.py %s -sample_name %s %s -out %s/%s/ " %( python_path , dnascan_dir , option_string , sample_name , input_file_string , out_dir , sample_name) )
        
        