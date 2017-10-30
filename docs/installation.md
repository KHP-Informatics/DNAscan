[TOC]

### Minimum requierements



- Ubuntu >= 14.04
- RAM: at least 4.2 Gb for the fast mode. 8 Gb for the other modes. 
- Space required by the installation: 
  - By the dependencies and DNscan: we reccomend at least 10Gb of free space for the whole DNAscan deplyment, including the dependencies. 
  - 8Gb for the reference human genome and the hisat2 index. An extra 5.5Gb for the bwa index are required if you want to use         the normal or intensive modes   
  - An extra 400Gb would be required if you wish to perform the annotation stage. This uses annovar databases which need to be dowloaded, such as CADD, Clinvar, dbSNP, etc
- Scratch space for usage: If you are performing the alignment stage and your input data is in fastq.gz format we recommend at 3 times the sice of your input data. E.g. your fastq.gz files are 100Gb, thus you would need 300Gb of free space. If you dont wish to performe alignment a proportion data-to-analyse:free-space of 1:1 would be enough. E.g. Input data is a 50Gb bam file, you would need only 50Gb free space on drive. 

### Obtaining

**Version:** 1.0

Please make sure all dependencies are installed before running DNAscan. Instrutions about how to install all dependencies are available in the following chapter. However a bash script to set up all dependencies is available in scripts.

To download DNAscan please use git to download the most recent development tree:

```bash
git clone https://github.com/snewhouse/DNA-NGS_scan.git
```

Once you have downloaded DNAscan, you can set up all available dependencies running the install_dependencies.sh script available in DNA-NGS_scan/scripts. This script will install all softaware dependencies as well as hg19 reference genome and its hisat2 and bwa indexes (these jobs runs in background and will finish after the scripts ends) as well as update paths.py. If you want your results being anotated you need to download Annovar. Since downloading Annovar requires registration, please register before using DNAscan [link](http://www.openbioinformatics.org/annovar/annovar_download_form.php) 

```bash
bash /path/to/DNAscan/scripts/install_dependencies.sh /path/to/set_up/directory /path/to/DNAscan/directory 

source ~/.bashrc

```

You can easly set up your running deployment od DNAscan using docker. For instructutions about how to install docker see following sections.

After installing docker run an Ubuntu image:

```bash
docker run -v /path/to/your/data_folder:/container/path/where/you/want/your/data [-p 8080:8080] -it [--storage-opt size=500G] ubuntu /bin/bash 

```

The -v option adds your data folder (assuming you have some data to run DNAscan on), -p mirrows the container port 8080 to your host port 8080. This will be necessary if you want to use the iobio services.
The --storage-opt size=NG option defines the maximum size of the container, setting it to N gigabytes. This number would depend on you plans. If you want to perform annotation, the databases used by DNAscan (clinvar,CADD,etc) are about 350G. We recomend N = 500 if you want to install the whole pipeline (including annotation). To this number you should add what needed for your analysis, e.g. if you are planning to download data, the size of you data to analyse etc. A workaround to this is to use the mirrowed host folder as outdir for your analysis. This folder does not have such size limit.  

IMPORTANT: to detach from the container without stopping it use Ctrl+p, Ctrl+q
IMPORTANT: when running DNAscan inside a docker container, if you want to use the iobio services (and for example upload your results into the gene.iobio platform), these would not be visible by your browser unless in the folder where they are is mirrowed on the host system. Considering the previous command to run an ubuntu image, the easiest way to do this would be using the folder where you imported the data inside the container (container/path/where/you/want/your/data) as out dir when running DNAsca. In this way the DNAscan results can could be found in /path/to/your/data_folder on the host system.

Then install git and download this repository and run the install_dependencies.sh script:

```bash
apt-get update

apt-get install git

git clone https://github.com/snewhouse/DNA-NGS_scan.git

cd DNA-NGS_scan

#if you are not interested in performing annotation, please # the annovar lines in the install_dependencies.sh script. This would avoid the download of about 350G od databses
bash scripts/install_dependencies.sh /path/to/set_up/directory /path/to/DNAscan/directory 

source ~/.bashrc

```

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



### Dependencies

#### List of dependencies

Fast mode pipeline (ideal if focusing on SNVs):

- Samtools 1.3
- HISAT2 2.1.0
- Freebayes 1.0.2
- Python 3.5
- Vcftools 0.1.13 
- Bedtools2 2.25
- Manta 1.2.0 (optional, needed only if interested in structural variants)
- ExpansionHunter 2.0.9 (optional, needed only if interested in known motif expansions)
- Bcftools 1.5 (optional, needed only if interested in performing custome variant filtering)
- Annovar "Version: $Date: 2016-02-01 00:11:18 -0800 (Mon, 1 Feb 2016)" (optional, needed only if interested in performing variant annotation)

Normal mode pipeline (better performance on indels and SVs):

- BWA 0.7.15 

Intensive mode pipeline (top performance on indels):

- Genome Analysis Toolkit 3.5 

Tools needed for generating graphical reports 

- RTG Tools 3.6.2 
- Multiqc 1.2 

Tools needed to allow an on-the-fly result interpretation 

- Gene.IoBio platform 2.1 
- Vcf.IoBio platform 
- Bam.IoBio platform

Tools needed for a container based deplyment 

- Docker 1.7.1
- Docker-compose 1.4.2
- Singularity 2.2 

#### How to install dependencies

##### Using Bioconda

For a fast and easy deployment of most dependencies we recomend the use of the Miniconda2 package.
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

Bioconda is a repository of binary bioinformatics tools which makes it very easy to install many open source software. 
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

Gene.iobio be found at the Tony Di Sera (https://github.com/tonydisera) github repo. Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
git clone https://github.com/tonydisera/gene.iobio.git
```

###### Vcf.iobio

Bam.iobio be found at the Tony Di Sera (https://github.com/tonydisera) github repo. Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
git clone https://github.com/tonydisera/vcf.iobio.io.git
```

###### Bam.iobio

Bam.iobio be found at the Chase Miller (https://github.com/chmille4) github repo. Please use git:

```bash
mkdir /path/to/your/iobio/
cd /path/to/your/iobio/
https://github.com/chmille4/bam.iobio.io.git
```

##### Download Annovar

A more complete documentation about how to set up and use annovar can be found [here](http://annovar.openbioinformatics.org/en/latest/). However, in the following, some brief instructions about how to get annovar and the nececcary databases to use DNAscan annotation.

The latest version of ANNOVAR (2017Jul16) can be downloaded [here](http://www.openbioinformatics.org/annovar/annovar_download_form.php) (registration required).

After you have downloaded Annovar:

```bash
tar xvfz annovar.latest.tar.gz
```

Now let's download the data bases needed for the DNAscan annotation step (assuming you want to work with hg19):

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

### Docker/Singularity

**Docker** is an open-source project that automates the deployment of applications inside software containers. Using containers to deploy our system and creating our analysis environment would allow us to make our work independent by the machine we work on. This would improve the reproducibility of our science, the portability and reliability of our deployments and avoid any machine specific issues. For this reason working using containers isn't just recommended but also makes things easier. Since docker is widely used and maintained we recommend it as container technology to use if possible. Unfortunately Docker does require sudo privileges to run its containers making its use difficult on HPC facilities.

**Singularity** is also a container project similar to Docker and does not require sudo privileges to run. This can be very important if you decide to use our framework on a machine for which you do not have such privileges. E.g. your institution HPC cluster. In this case you create your docker deployment locally and then converting the docker image into a singularoty image using this [script](https://github.com/KHP-Informatics/MNDA-DataManagement-System/tree/master/docker2singularity)

```bash 
$ ./docker2singularity.sh  [-m \"/mount_point1 /mount_point2\"] docker_image_name
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

- test <br />
  python run_analysis.py -format input_file_fomat -reference hg19 -in input_file -out out_dir -SV -BED<br />



### How to download the reference genome

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

### How to download the NCBI microbes DBs

#### Virus

Copy and paste in your command line the following commands to download the whole NCBI database of complete viral genome

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

Copy and paste in your command line the following commands to download the whole NCBI database of complete bacterial genome (this might take a long time)

```bash 
mkdir /path/to/wherever/bacteria_db

cd /path/to/wherever/bacteria_db

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.*.1.genomic.fna.gz

for i in $(ls | grep -e genomic.fna.gz -e bacteria); do gzip -d $i ; rm $i ; done 

for i in $(ls | grep -e genomic.fna -e bacteria); do cat $i >> bacteria_db.fa ; rm $i ; done  
```

#### Fungi

Copy and paste in your command line the following commands to download the whole NCBI database of complete bacterial genome (this might take a long time)

```bash 
mkdir /path/to/wherever/fungi_db

cd /path/to/wherever/fungi_db

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/fungi.*.1.genomic.fna.gz

for i in $(ls | grep -e genomic.fna.gz -e fungi); do gzip -d $i ; done 

for i in $(ls | grep -e genomic.fna -e fungi); do cat $i >> fungi_db.fa ; rm $i ; done  
```

### How to index the reference genome (or a microbes DB)

#### HISAT2

```bash
./path/to/hisat/hisat2-build /path/to/reference/file/reference_genome.fa index_base
```

E.g. If the reference genome is the file hg19.fa, located in /home/dataset/ and the hisat2-build binary is located in /home/bin/, the command line would be: 

```bash
./home/bin/hisat2-build /home/dataset/hg19.fa hg19
```



#### BWA

```bash
./path/to/bwa/bwa index /path/to/reference/file/reference_genome.fa
```

E.g. If the reference genome is the file hg19.fa, located in /home/dataset/ and the bwa binary is located in /home/bin/, the command line would be: 

```bash
./home/bin/bwa index /home/dataset/hg19.fa
```

