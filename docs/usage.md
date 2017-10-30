[toc]

### Usage

IMPORTANT: DNAscan.py is the main script performing the analyses. It must be in the same folder as paths.py. Before running DNAscan please modify paths.py to match your dependencies deplyment.
IMPORTANT2: All paths in DNAscan end with "/"

Its basic use requires the following options:

```bash
  -format FORMAT        options are bam, sam, fastq, vcf [string] 
  -reference REFERENCE  options are hg19, hg38 [string]
  -in INPUT_FILE        input file [string]
  -out OUT              path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/

```

 The desired pipeline stages are performed according to the optional arguments selected:

```bash
  -filter_string FILTER_STRING  bcftools filter string, eg GQ>20 & DP>10 (Default = "")
  -iobio                if this flag is set iobio services will be started at the end of the analysis (Default = "False")
  -alignment            if this flag is set the alignment stage will be performed (Default = "False")
  -expansion            if this flag is set DNAscan will look for the expansions described in the jason folder described in paths.py  (Default = "False") 
  -SV                   if this flag is set the structural variant calling stage will be performed (Default = "False") 
  -virus                if this flag is set DNAscan will perform viral scanning (Default = "False")  
  -bacteria             if this flag is set DNAscan will perform bacteria scanning (Default = "False") 
  -custom_microbes      if this flag is set DNAscan will perform a customized microbe scanning according to the provided microbe data base in paths.py (Default = "False")  
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

Fast mode uses Hisat2 and Freebayes to quickly align and call variants. It is ideal if you are focusing your analysis on single nucleotyde variants. Normal mode performs an alignment refinement using BWA on selected reads. This step improves the alignment of soft-clipped reads and reads containing small indels. It is suggested if your focus is on structural variants. Intensive mode adds to the pipeline a further indel calling using GATK Haplotype Caller which improves the performance on small indels. If your analysis focus on the discovery of non human material (e.g. viruses or bacteria) in your sequencing data, fast mode is recomended for known microbes discovery while normal mode improves the discovery or divergent microbes. A detailed description of the 3 modes can be found in the [DNAscan paper](link).  

Finally, a set of optional arguments can be used to customise the analysis:

```bash
-RG RG                if this flag is set the alignment stage will add the provided in paths.py read group (Default = "False")
-paired PAIRED        options are 1 for paired end reads and 0 for single end reads (Default = "1")
-vcf VCF              complementary vcf file
-in2 INPUT_FILE2      input file 2, for paired end reads only (usually fastq file)
-BED                  restrict the analysis to the regions in the bed file (Default = "False") 
-sample_name SAMPLE_NAME  specify sample name [string] (default = "sample")
-debug                if this flag is set DNAscan will not delete intermediete and temporary files (Default = "False")
```

#### Usage example

Let's assume we have human paired end whole exome sequening data in two fastq files and want to perform snvs/indels calling vs hg19, annotation and explore the results using the iobio services. The DNAscan command line would be:

```bash
python3 /path/to/DNAscan/scripts/DNAscan.py -format fastq -in data1.fq.gz -in2  data2.fq.gz -reference hg19 -alignment -variantcalling -annotation -iobio -out /path/to/outdir -mode fast
```

Using the sequencing data provided in the data folder:

```bash
 cd /path/to/DNAscan_main_dir
 
python3 scripts/DNAscan.py -format fastq -in data/test_data.1.fq.gz -in2 data/test_data.2.fq.gz -reference hg19 -alignment-variantcalling -annotation -iobio - out outdir/ -mode fast -BED
```

IMPORTANT: All paths in DNAscan end with "/"