![alt text](https://github.com/KHP-Informatics/DNAscan/blob/master/DNAscan_logo.001.jpeg)

# DNAscan
```diff
+TO BE NOTED: We are always working to improve DNAscan so any bug reports, suggestions, feedbacks would be highly welcome. 
```
## Table of Contents
1. [Introduction](#introduction)
2. [Citation](#citation)
3. [Documentation](#documentation)
	* [Minimum requirements](#minimum-requirements)
	* [Obtaining](#obtaining)
	* [Usage](#usage)
	* [ALSgeneScanner](#alsgenescanner)
	* [Output](#output)
	* [How to download the reference genome](#how-to-download-the-reference-genome)
	* [How to download the NCBI microbe DBs](#how-to-download-the-ncbi-microbe-dbs)
	* [How to index the reference genome or a microbe DB](#how-to-index-the-reference-genome-or-a-microbe-db)
	* [Graphic user interface](#Graphic-user-interface)
	* [Dependencies](#dependencies)
	* [Docker and Singularity](#docker-and-singularity)	
4. [Core Contributors](#core-contributors)
5. [Contributing](#contributing)
6. [Licence](#licence)


## Introduction

DNAscan is a fast and efficient bioinformatics pipeline that allows for the analysis of DNA Next Generation sequencing data, requiring very little computational effort and memory usage. DNAscan can analyse 40x whole genome NGS data in ~8 hours, using as little as 8 cpus and 16 Gbs of RAM while guaranteeing a very high performance. We do this by exploiting state-of-the-art bioinformatics tools. DNAscan can screen your DNA NGS data for single nucleotide variants, small indels, structural variants, repeat expansions, viral (or any other organism’s) genetic material. Its results are annotated using a wide range of databases including ClinVar, EXAC, dbSNP and CADD and uploaded onto the gene.iobio platform for an on-the-fly analysis/interpretation.

![alt text](https://github.com/KHP-Informatics/DNAscan/raw/master/DNAscan_pipeline_pic.jp2)

Figure 1. Central panel: Pipeline overview. DNAscan accepts sequencing data, and optionally variant files. The pipeline firstly performs an alignment step (details in the left panel), followed by a customizable data analysis protocol (details in the right panel). Finally, results are annotated and a user friendly report is generated. Right panel: detailed description of the post alignment analysis pipeline (intensive mode). Aligned reads are used by the variant calling pipeline (Freebayes + GATK HC + Bcftools), both aligned and unaligned reads are used by Manta and ExpensionHunter (for which known repeat description files have to be provided)  to look for structural variants. The unaligned reads are mapped to a database of known viral genomes (NCBI database) to screen for the presence of their genetic material in the input sequencing data. Left panel: Alignment stage description. Raw reads are aligned with HISAT2. Resulting soft/hard-clipped reads and unaligned reads are realigned with BWA mem and then merged with the others using samtools.


## Citation

[Alfredo Iacoangeli et al. DNAscan: a fast, computationally and memory efficient bioinformatics pipeline for analysis of DNA next-generation-sequencing data. bioRxiv, 2018](https://www.biorxiv.org/content/biorxiv/early/2018/02/18/267195.full.pdf)

## Documentation

### Minimum requirements

- Ubuntu >= 14.04.
- RAM: if processing WGS data, at least 14 Gb for the fast mode. 16 Gb for the other modes. 
- Space required by the installation: 
  - By the dependencies and DNAscan: we recommend at least 10Gb of free space for the whole DNAscan deployment, including the dependencies. 
  - 8Gb for the reference human genome and the hisat2 index. An extra 5.5Gb for the bwa index are required if you want to use         the normal or intensive modes.   
  - The annotation step makes use of Annovar and selected databases, the size of which can range from tens of megabytes to hundreds of gigabytes (e.g. CADD database). If the user wishes to perform the annotation step they must take this into account.
- Scratch space for usage: If you are performing the alignment stage and your input data is in fastq.gz format we recommend using at least 3 times the size of your input data. E.g. your fastq.gz files are 100Gb, thus you would need 300Gb of free space. If you don't wish to perform the alignment, a proportion of data-to-analyse:free-space of 1:1 would be enough. E.g. Input data is a 50Gb bam file, you would need only 50Gb of free space. 

### Obtaining

**Version:** 0.1

Please make sure all dependencies are installed before running DNAscan. Instructions on how to install all dependencies are available in the following chapter. However a bash script to set up all dependencies (Annovar and GATK need a manual registration and download step) is available in scripts.

#### Local Deployment

To obtain DNAscan please use git to download the most recent development tree:

```bash
git clone https://github.com/KHP-Informatics/DNAscan.git
```

Once you have downloaded DNAscan, you can set up all needed dependencies running the install_dependencies.sh script available in DNAscan/scripts. Before running install_dependencies.sh please download and uncompress Annovar by registering at the following [link](http://download.openbioinformatics.org/annovar_download_form.php) and download GATK 3.8 at the following [link](https://github.com/broadgsa/gatk/releases). Install_dependencies.sh will install all software dependencies as well as hg19 reference genome and its hisat2 and bwa indexes (these jobs run in the background and will finish after the script ends) as well as update paths_and_configs.py. 

```bash

bash scripts/install_dependencies.sh /path/to/set_up/directory/ /path/to/DNAscan/directory/ /path/to/annovar/directory/ /path/to/GATK3.8/download.tar.bz2 $num_threads

source ~/.bashrc

```

#### Obtain with Docker

IMPORTANT: if you want to use the intensive mode or the annotation step of DNAscan you need to register and download both Annovar and GATK 3.8. Remember to mirror into the container the folder where you have deployed these using the -v flag or copy them inside the container using the docker cp command as described below.
The easiest way to get started with DNAscan is to use its Docker image:

```bash
sudo docker run  -v /path/to/your/data_folder:/container/path/where/you/want/your/data [-p 8080:8080] --storage-opt size=50G  -it compbio/dnascan /bin/bash
```
Please set the needed container size taking into account the "Minimum requirements".

The -v option adds your data folder (assuming you have some data to run DNAscan on), -p mirrors the container port 8080 to your host port 8080. This will be necessary if you want to use the iobio services.
The --storage-opt size=Ngigas option defines the maximum size of the container, setting it to N gigabytes. This number would depend on you plans. If you want to perform annotation, the databases used by DNAscan (clinvar,CADD,etc) are about 350G. We recommend N = 500 if you want to install the whole pipeline (including annotation). To this number you should add what you need for your analysis, e.g. if you are planning to download data, the size of your data to analyse etc. A way to workaround this is to use the mirrored host folder as outdir for your analysis and as annovar folder. This folder does not have a size limit. 

IMPORTANT: To detach from the container without stopping it, use Ctrl+p, Ctrl+q.
IMPORTANT: When running DNAscan inside a docker container, if you want to use the iobio services (and for example upload your results into the gene.iobio platform), these would not be visible by your browser unless they are in the folder which is mirrored on the host system. Considering the previous command to run an ubuntu image, the easiest way to do this would be to use the folder where you imported the data inside the container (/container/path/where/you/want/your/data) as an outdir when running DNAsca. In this way the DNAscan results can be found in /path/to/your/data_folder on the host system.

If you want to add data to the container while this is already running you can use the docker cp command. First detach from the container without stopping it using Ctrl+p, Ctrl+q, then cp your data inside the container:

```bash
docker cp [OPTIONS] /path/to/your/data CONTAINER_ID:/container/path/where/you/want/your/data
```
and execute a bash shell inside your container:

```bash
docker exec -it CONTAINER_ID /bin/bash
```
The container ID can be found using the ps command:

```bash
docker ps 
```


#### Building a Docker container from scratch

You can easily set up your running deployment of DNAscan using docker. For instructutions about how to install docker see following.

After installing docker run an Ubuntu image:


```bash

docker run -v /path/to/your/data_folder:/container/path/where/you/want/your/data [-p 8080:8080] -it [--storage-opt size=500G] ubuntu /bin/bash 

```

Then install git, download this repository and run the install_dependencies.sh script:

```bash

apt-get update

apt-get install git

git clone https://github.com/KHP-Informatics/DNAscan.git

cd DNAscan

#By dafault install_dependencies.sh downloads the following Annovar databases: Exac, Refgene, Dbnsfp, Clinvar and Avsnp. If you wish to download the CADD database (about 350G) please uncomment the appropiete line. If you are not interested in performing annotation, please # the annovar lines in the install_dependencies.sh script. 

bash scripts/install_dependencies.sh /path/to/set_up/directory/ /path/to/DNAscan/directory/ /path/to/annovar/directory/ /path/to/GATK3.8/download.tar.bz2 $num_threads

source ~/.bashrc

```

### Usage

IMPORTANT: DNAscan.py is the main script performing the analyses. It must be in the same folder as paths_and_configs.py. Before running DNAscan please modify paths_and_configs.py to match your dependencies deplyment.
IMPORTANT2: All paths in DNAscan end with "/"

Its basic use requires the following options:

```bash

  -format FORMAT        options are bam, sam, fastq, vcf [string] 
  -reference REFERENCE  options are hg19, hg38, grch37 and grch38 [string]
  -in INPUT_FILE        input file [string]
  -in2 INPUT_FILE       second input file (for paired end reads in fastq format only) [string]
  -out OUT              path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/

 ```
 The desired pipeline stages are performed according to the optional arguments selected:
 
 ```bash
  -filter_string FILTER_STRING  bcftools filter string, eg GQ>20 & DP>10 (Default = "")
  -iobio                if this flag is set iobio services will be started at the end of the analysis (Default = "False")
  -alignment            if this flag is set the alignment stage will be performed (Default = "False")
  -expansion            if this flag is set DNAscan will look for the expansions described in the json folder described in paths_and_configs.py  (Default = "False"). It requires a path to a folder containing the json repeat-specification files to be specified in paths_and_configs.py
  -SV                   if this flag is set the structural variant calling stage will be performed (Default = "False") 
  -virus                if this flag is set DNAscan will perform viral scanning (Default = "False")  
  -bacteria             if this flag is set DNAscan will perform bacteria scanning (Default = "False") 
  -custom_microbes      if this flag is set DNAscan will perform a customized microbe scanning according to the provided microbe data base in paths_and_configs.py (Default = "False")  
  -variantcalling       if this flag is set DNAscan will perform snv and indel calling (Default = "False")  
  -annotation           if this flag is set DNAscan will annotate the found variants (Default = "False")  
  -results_report       if this flag is set DNAscan generate a results report (Default = "False")  
  -alignment_report     if this flag is set DNAscan generate an alignment report (Default = "False") 
  -sequencing_report    if this flag is set DNAscan generate a report describing the input sequencing data (Default = "False") 
  -calls_report         if this flag is set DNAscan generate a report describing the found snvs and indels (Default = "False")
  -rm_dup               if this flag is set DNAscan will remove duplicates from the sample bam file (Default = "False")

```

Also, one of the three analysis modes can be chosen with the -mode option:

 ```bash
-mode  MODE            options are fast, normal, intensive [string] (default = "fast")

```
Fast mode uses Hisat2 and Freebayes to quickly align and call variants. It is ideal if you are focusing your analysis on single nucleotyde variants. Normal mode performs an alignment refinement using BWA on selected reads. This step improves the alignment of soft-clipped reads and reads containing small indels. It is recommended if your focus is on structural variants. Intensive mode adds a further indel calling step to the pipeline using GATK Haplotype Caller which improves the performance on small indels. If your analysis focuses on the discovery of non human material (e.g. viruses or bacteria) in your sequencing data, please note that that the selceted mode does not affect this step. A detailed description of the 3 modes can be found in the [DNAscan paper](https://www.biorxiv.org/content/biorxiv/early/2018/02/18/267195.full.pdf).  

Finally, a set of optional arguments can be used to customise the analysis:

 ```bash
-RG RG                if this flag is set the alignment stage will add the provided in paths_and_configs.py read group (Default = "False")
-paired PAIRED        options are 1 for paired end reads and 0 for single end reads (Default = "1")
-vcf VCF              complementary vcf file
-in2 INPUT_FILE2      input file 2, for paired end reads only (usually fastq file)
-BED                  restrict the analysis to the regions in the bed file (Default = "False") 
-sample_name SAMPLE_NAME  specify sample name [string] (default = "sample")
-debug                if this flag is set DNAscan will not delete intermediete and temporary files (Default = "False")
```


#### How to make the read mappers and variant callers use custom options

The file paths_and_configs.py in DNAscan_installation_folder/scripts can be used for this purpose. In paths_and_configs.py there is a section where a string of options can be passed to the read mappers (Hisat2 and BWA) and the variant callers (Freebayes and GATK HC):

 ```bash
 #custom tool options

hisat_custom_options = ""

bwa_custom_options = ""

freebayes_custom_options = ""

gatk_HC_custom_options = ""
```

For example, if we want the callers to discard reads with mapping quality < 15, we can modify as follows:


```bash
freebayes_custom_options = "--min-mapping-quality 15"

gatk_HC_custom_options = "-mmq 15"
```

For a complete list of options that can be used with these tools please see the following links:

[Hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
[BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
[Freebayes](https://github.com/ekg/freebayes)
[GATK HC](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/?v=3)

#### Usage example

Let's assume we have human paired end whole exome sequening data in two fastq files and want to perform snvs/indels calling vs hg19, annotation and explore the results using the iobio services. The DNAscan command line would be:

 ```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out /path/to/outdir/ -mode fast
```
Using the sequencing data provided in the data folder:

 ```bash
cd /path/to/DNAscan_main_dir
 
python3 scripts/DNAscan.py -format fastq -in data/test_data.1.fq.gz -in2 data/test_data.2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out outdir/ -mode fast -BED
```
IMPORTANT: All paths in DNAscan end with "/"

#### Applying filters to the called variants

DNAscan uses bcftools to filter the variants. It only selects variants for which the expression provided is true. A few examples follow. Complete details about how to write an expression can be found in the bcftool manual [LINK](https://samtools.github.io/bcftools/bcftools.html#expressions).  

Filtering calls for which genotype quality is > then 30:

```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out /path/to/outdir/ -mode fast -filter "GQ>30"
```
Filtering calls for which genotype quality is > then 30 and depth >5:

```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out /path/to/outdir/ -mode fast -filter 'GQ>30 && DP>5'
```


#### Looking for repeat expansions

Let's assume we have human paired end whole exome sequening data in two fastq files and want to perform snvs/indels calling vs hg19 and scan for some specific repeat expansions. The DNAscan command line would be:

 ```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -expansion -out /path/to/outdir/ -mode fast
```
Note that the json repeat-specification files must be specified in paths_and_configs.py. For a guide on how to create a json repeat-specification file see the this [LINK](https://github.com/Illumina/ExpansionHunter/wiki/Inputs) to the ExpantionHunter instructions. Two examples of such files. for the C9orf72 repeat and ataxin2 repeat are in DNAscan/repeats folder.

#### Running DNAscan on a subregion of the human region

##### Whole exome

You can restrict the analysis to the whole exome only. You just need to add the -exome flag to the command line:

 ```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -out /path/to/outdir/ -exome
```

##### Any genome regions

You can restrict the analysis to any specific region of the genome by providing DNAscan with a bed file. In this case you need to use the -BED flag when running the pipeline and specify the path to the bed file in paths_and_configs.py

 ```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alignment -variantcalling -out /path/to/outdir/ -BED
```
##### List of genes

You can restrict the analysis to a list of genes. In this case you only need to generate a list of genes (see example in DNAscan/data) and specify the path to it in paths_and_configs.py. No flags need to be used in this case and DNAscan will read (non-case-sensitive) the list and restrict the analsyis to the appropriate genome regions. 

 ```bash
cat list_of_genes.txt

fus
PFN1
SARM1
SCFD1
sod1
```

#### Running DNAscan on a list of BAMs/FASTQs

Let's assume we have human paired end whole exome sequening data in two fastq files per sample and want to perform snvs/indels calling vs hg19. In this case we would use the analyse_list_of_samples.py script to run DNAscan:

 ```bash
 cd /path/to/DNAscan_main_dir
 
python3 scripts/analyse_list_of_samples.py -option_string "-format fastq -reference hg19 -alignment -variantcalling -mode fast" -out_dir /path/to/where/you/want/the/analysis/to/take/place/ -sample_list extras/fastq_sample_list.txt -format fastq -paired 1 
```
The analyse_list_of_samples.py takes as input a file containing the input file of one sample per line (tab separated in case of 2 paired end reads fastq files per sample), the option you want to pass to DNAscan between quotation marks, the input file format and whether or not (1 or 0) the data is paired-end reads. 2 sample lists of input files are in DNAscan/extras

### ALSgeneScanner

ALSgeneScanner is a pipeline designed for the analysis of NGS data of ALS patients. It perfoms alignment, variant calling, structural variant callin and repeat expansion calling as well as variant annotation using Annovar. It restricts the analysis to a subset of genes (150) which have been shown to be associated with ALS. A complete list of the included genes is available as a [Google Spreadsheet](https://goo.gl/mngrMC). It also prioritize variants according to the scientific evidence of the gene association and the effect prediction of the variant. To run ALSgeneScanner the user has to use the appropriate flag as in the following example. At present this can only be done using the reference genome hg19.

Usage example:

```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2 data2.fq.gz -reference hg19 -alsgenescanner -out /path/to/outdir/ -mode fast
```

RefGene, ClinVar and dbnsfp30a are nessary to run ALSgeneScanner.

Please use the following commands to download the appropiate Annovar databases:

```bash
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a /path/to/annovar/database/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene /path/to/annovar/database/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20170905 /path/to/annovar/database/
```


### Output

DNAscan output tree:

```bash

./$out_dir-| # this is the folder given to DNAsca using the -out flag. It will contain the aligned sequecing data ($sample_name.bam) as well as some temporanery files
           |
           |-results # This will contain the output of the analyses. E.g. $sample_name_sorted.vcf.gz , $sample_name_SV.vcf.gz, virus_results.txt, etc  
           |
           |-reports # If any report flags is used, this folder will contain the reports. E.g. $sample_name_vcfstats.txt if the -calls_report flag is used
           |
           |-logs # Logs files are generated in this folder
 
```          
           
### How to download the reference genome

DNAscan can be used with the following human genome reference builds: hg19, hg38, grch37, grch38. The used build can be specified when running DNAscan by using the option -reference (options are hg19, hg38, grch37 and grch38). However, DNAscan uses Annovar to annotate the variants and Annovar is only compatible with hg19 and hg38 at present. Therefore, if the user wants to use grch37 and grch38 the pipeline will skip the annotation step. We are currently working on this and it will be possible to run DNAscan on grch37 and grch38 without limitations. 

#### hg38

```bash 
mkdir /path/to/wherever/hg38

cd /path/to/wherever/hg38

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gzip -d hg38.fa.gz

samtools faidx hg38.fa
```

#### hg19
```bash 
mkdir /path/to/wherever/hg19

cd /path/to/wherever/hg19

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

tar -zxvf chromFa.tar.gz

for i in chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrY.fa chrX.fa chrM.fa; do cat $i >> hg19.fa ; rm $i ; done

rm chr*

samtools faidx hg19.fa
```
#### grch37
```bash 
mkdir /path/to/wherever/grch37

cd /path/to/wherever/grch37

for i in {1,2,3,4,5,6,7,8,9}; do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/CHR_0$i/hs_ref_GRCh37.p5_chr$i.fa.gz; zcat hs_ref_GRCh37.p5_chr$i.fa.gz >> grch37.fa ; rm hs_ref_GRCh37.p5_chr$i.fa.gz; done


for i in {10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT}; do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/CHR_$i/hs_ref_GRCh37.p5_chr$i.fa.gz; zcat hs_ref_GRCh37.p5_chr$i.fa.gz >> grch37.fa ; rm hs_ref_GRCh37.p5_chr$i.fa.gz; done

samtools faidx grch37.fa

```
#### grch38
```bash 
mkdir /path/to/wherever/grch38

cd /path/to/wherever/grch38

for i in {1,2,3,4,5,6,7,8,9}; do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/CHR_0$i/hs_ref_GRCh38.p12_chr$i.fa.gz; zcat hs_ref_GRCh38.p12_chr$i.fa.gz >> grch38.fa ; rm hs_ref_GRCh38.p12_chr$i.fa.gz; done


for i in {10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT}; do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/CHR_$i/hs_ref_GRCh38.p12_chr$i.fa.gz; zcat hs_ref_GRCh38.p12_chr$i.fa.gz >> grch38.fa ; rm hs_ref_GRCh38.p12_chr$i.fa.gz; done

samtools faidx grch38.fa

```

### How to download the NCBI microbe DBs

#### Virus

Copy and paste the following commands in your command line to download the whole NCBI database of complete viral genomes.

```bash 
mkdir /path/to/wherever/virus_db

cd /path/to/wherever/virus_db

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz

gzip -d viral.1.1.genomic.fna.gz

gzip -d viral.2.1.genomic.fna.gz

cat viral.1.1.genomic.fna >> virus_db.fa

cat viral.2.1.genomic.fna >> virus_db.fa

rm viral.1.1.genomic.fna viral.2.1.genomic.fna 
```


#### Bacteria

Copy and paste the following commands in your command line to download the whole NCBI database of complete bacterial genomes (this might take a long time).

```bash 
mkdir /path/to/wherever/bacteria_db

cd /path/to/wherever/bacteria_db

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.*.1.genomic.fna.gz

for i in $(ls | grep -e genomic.fna.gz -e bacteria); do gzip -d $i ; rm $i ; done 

for i in $(ls | grep -e genomic.fna -e bacteria); do cat $i >> bacteria_db.fa ; rm $i ; done  
```

#### Fungi

Copy and paste the following commands in your command line to download the whole NCBI database of complete fungi genomes. 

```bash 
mkdir /path/to/wherever/fungi_db

cd /path/to/wherever/fungi_db

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/fungi.*.1.genomic.fna.gz

for i in $(ls | grep -e genomic.fna.gz -e fungi); do gzip -d $i ; done 

for i in $(ls | grep -e genomic.fna -e fungi); do cat $i >> fungi_db.fa ; rm $i ; done  
```

### How to index the reference genome or a microbe DB

#### HISAT2
```bash
/path/to/hisat/hisat2-build /path/to/reference/file/reference_genome.fa index_base
```
E.g. If the reference genome is the file hg19.fa, located in /home/dataset/ and the hisat2-build binary is located in /home/bin/, the command line would be: 

```bash
/home/bin/hisat2-build /home/dataset/hg19.fa hg19
```



#### BWA 
```bash
/path/to/bwa/bwa index /path/to/reference/file/reference_genome.fa
```

E.g. If the reference genome is the file hg19.fa, located in /home/dataset/ and the bwa binary is located in /home/bin/, the command line would be: 

```bash
/home/bin/bwa index /home/dataset/hg19.fa
```

### Graphic user interface

We made a graphic user interface (GUI) for DNAscan to enable a wider audience to use it. The GUI source code is in DNAscanUGUI.zip. We used [UGUI](http://ugui.io/) to develop our GUI. To get it working please follow a few simple steps described below, however, detailed documentation and instruction about how to get it started can be found [here](http://ugui.io/tutorials/getting-started.html). If you do not have Node.js installed please get it from [here](https://nwjs.io/).

```bash
export PATH=/path/to/DNAscan/scripts/:$PATH

unzip DNAscanUGUI.ip
cd UGUI
npm install
npm start
```


### Dependencies

#### List of dependencies

Fast mode pipeline (ideal if focusing on SNVs):
* Samtools >= 1.5
* HISAT2 >= 2.1.0
* Freebayes >= 1.0.2
* Python >= 3
* Vcftools >= 0.1.13 
* Bedtools2 >= 2.25
* Samblaster >= 0.1.24
* Sambamba >= 0.6.6
* Manta 1.2.0 (optional, needed only if interested in structural variants)
* ExpansionHunter >= 2.0.9 (optional, needed only if interested in known motif expansions)
* Bcftools >= 1.3 (optional, needed only if interested in performing custome variant filtering and calls report)
* Annovar "Version >= $Date: 2016-02-01 00:11:18 -0800 (Mon, 1 Feb 2016)" (optional, needed only if interested in performing variant annotation)

Normal mode pipeline (better performance on indels and SVs):
* BWA 0.7.15 

Intensive mode pipeline (top performance on indels):
* Genome Analysis Toolkit 3.8 

Tools needed for generating graphical reports 
* RTG Tools >= 3.6.2 
* Multiqc >= 1.2 

Tools needed to allow an on-the-fly result interpretation 
* Gene.IoBio platform >= 2.1 
* Vcf.IoBio platform 
* Bam.IoBio platform

Tools needed for a container based deplyment 
* Docker >= 1.7.1
* Singularity >= 2.2 

#### How to install dependencies

##### Using Bioconda

For a fast and easy deployment of most dependencies we recommend the use of the Miniconda2 package.
To download and install the latest Miniconda2 package, which contains the conda package manager:

```bash
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh -b -p /path/to/Miniconda2/installation/directory
```

In the file /path_to_your_home_dir/.bashrc add the following line:

```bash
export PATH=/path/to/Miniconda2/installation/directory/Miniconda2/bin:$PATH 
```
Bioconda is a repository of binary bioinformatics tools which makes it very easy to install many open source softwares. 
To install the needed dependencies with Bioconda:

```bash
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

And then, if you want to install samtools for example:

```bash
conda install samtools
```

##### Download iobio services

###### Gene.iobio

Gene.iobio can be found at the Tony Di Sera github repo (https://github.com/tonydisera). Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
git clone https://github.com/tonydisera/gene.iobio.git
```

###### Vcf.iobio

Bam.iobio can be found at the Tony Di Sera (https://github.com/tonydisera) github repo. Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
git clone https://github.com/tonydisera/vcf.iobio.io.git
```

###### Bam.iobio

Bam.iobio can be found at the Chase Miller (https://github.com/chmille4) github repo. Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
https://github.com/chmille4/bam.iobio.io.git
```

##### Download Annovar

A more complete documentation about how to set up and use Annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/). However, in the following, there are some brief instructions on how to get Annovar and the nececcary databases to use DNAscan annotation.

The latest version of Annovar can be downloaded [here](http://www.openbioinformatics.org/annovar/annovar_download_form.php) (registration required).

After you have downloaded Annovar:

```bash
tar xvfz annovar.latest.tar.gz
```
Now let's download some data bases for the DNAscan annotation step (assuming you want to work with hg19):

```bash
cd annovar

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20170130 humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd humandb/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd humandb/

```

Finally, we have to update paths_and_configs.py appropriately. Details about how to download and use the Annovar databases can be found [here](http://annovar.openbioinformatics.org/en/latest/). For the above databases:

```bash
annovar_protocols = "refGene,dbnsfp30a,clinvar_20170130,avsnp147,cadd"

annovar_operations = "g,f,f,f,f"

```

### Docker and Singularity 

**Docker** is an open-source project that automates the deployment of applications inside software containers. Using containers to deploy our system and creating our analysis environment would allow us to make our work independent of the machine we work on. This would improve the reproducibility of our science, the portability and reliability of our deployments and avoid any machine specific issues. For this reason working using containers isn't just recommended but also makes things easier. Since docker is widely used and maintained we recommend it as container technology to use if possible. Unfortunately Docker does require sudo privileges to run its containers making its use difficult on HPC facilities.

**Singularity** is also a container project similar to Docker and does not require sudo privileges to run. This can be very important if you decide to use our framework on a machine for which you do not have such privileges. E.g. your institution HPC cluster. In this case you can use Singularity to convert the docker image into aSsingularity image and run a bash shell in the resulting Singularity container:

```bash 
$ singularity shell docker://compbio/dnascan
```
After starting the bash shell inside the singularity container you can find a working deployment of DNAscan in /DNAscan

```bash 
$ Singularity.dnascan> cd /DNAscan
$ Singularity.dnascan> cat /DNAscan/docker/welcome_message.txt
```


##### Docker setup

###### ubuntu

[LINK](https://store.docker.com/editions/community/docker-ce-server-ubuntu)

###### MAC

[LINK](https://store.docker.com/editions/community/docker-ce-desktop-mac)

###### Windows

[LINK](https://store.docker.com/editions/community/docker-ce-desktop-windows)

##### Singularity setup

###### Linux

[LINK](http://singularity.lbl.gov/install-linux)

###### MAC

[LINK](http://singularity.lbl.gov/install-mac)


## Core Contributors
- [Dr Alfredo Iacoangeli](alfredo.iacoangeli@kcl.ac.uk), UK
- [Dr Stephen J Newhouse](stephen.j.newhouse@gmail.com), UK

### Contributors
- [NAME](email), [COUNTRY OR AFFILIATION]

For a full list of contributors see [LINK](./CONTRIBUTORS.md)

## Contributing

Here’s how we suggest you go about proposing a change to this project:

1. [Fork this project][fork] to your account.
2. [Create a branch][branch] for the change you intend to make.
3. Make your changes to your fork.
4. [Send a pull request][pr] from your fork’s branch to our `master` branch.

Using the web-based interface to make changes is fine too, and will help you
by automatically forking the project and prompting to send a pull request too.

[fork]: https://help.github.com/articles/fork-a-repo/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[pr]: https://help.github.com/articles/using-pull-requests/


## Licence 
- [MIT](./LICENSE.txt)


*********



<p align="center">
  Core Developers funded as part of:</br> 
  <a href="https://www.mndassociation.org/">MNDA</a></br> 
  <a href="http://www.maudsleybrc.nihr.ac.uk/">NIHR Maudsley Biomedical Research Centre (BRC), King's College London</a></br>
  <a href="http://www.ucl.ac.uk/health-informatics/">Farr Institute of Health Informatics Research, UCL Institute of Health Informatics, University College London</a>
</p>


![alt text](https://raw.githubusercontent.com/KHP-Informatics/MND-DataManagementAnalysis-System/master/funder_logos.001.jpeg)
