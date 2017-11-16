If you are reading this file you have already installed all DNAscan dependencies except GATK 3.8 and Annovar.
These two software require registration to be used. Please register and download Annovar at http://www.openbioinformatics.org/annovar/annovar_download_form.php and GATK 3.8 at https://software.broadinstitute.org/gatk/download/ .

After succefully registering and downloading the two software you can either change the appropiate paths in paths.py and download the needed Annovar databases (see below) or run the following command:

bash /path/to/DNAscan/scripts/install_annovar_and_gatk.sh /path/to/annovar/ /path/to/gatk.bz2/dir/  

To download the Annovar databases use the following commands replacing the appropiate paths:

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd $DNAscan_INSTALL_DIR/humandb/ &

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20170130 $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 $DNAscan_INSTALL_DIR/humandb/

