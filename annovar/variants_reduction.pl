#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Copy;
use File::Basename;

our $REVISION = '$Revision: 871849c7dd7b3842a6f3acecc3ce47beb5e1e920 $';
our $DATE =	'$Date: 2015-12-14 13:51:19 -0800 (Mon, 14 Dec 2015) $';  
our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org> $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $remove, $checkfile, $dispensable, $genetype, $aaf_threshold, $maf_threshold, $protocol, $operation, $genericdbfile,
	$ljb_sift_threshold, $ljb_pp2_threshold, $ljb2_sift_threshold, $ljb2_pp2hvar_threshold, $ljb2_pp2hdiv_threshold, $argument);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'remove'=>\$remove,
	'checkfile!'=>\$checkfile, 'dispensable=s'=>\$dispensable, 'genetype=s'=>\$genetype, 'maf_threshold=f'=>\$maf_threshold,
	'aaf_threshold=f'=>\$aaf_threshold, 'protocol=s'=>\$protocol, 'operation=s'=>\$operation, 'genericdbfile=s'=>\$genericdbfile, 
	'ljb_sift_threshold=f'=>\$ljb_sift_threshold, 'ljb_pp2_threshold=f'=>\$ljb_pp2_threshold, 'ljb2_sift_threshold=f'=>\$ljb_sift_threshold, 'ljb2_pp2hvar_threshold=f'=>\$ljb_pp2_threshold, 'ljb2_pp2hdiv_threshold=f'=>\$ljb2_pp2hdiv_threshold,
	'argument=s'=>\$argument) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");


($queryfile, $dbloc) = @ARGV;

$^O eq 'MSWin32' and die "Error: this program does not work in Microsoft Windows operating system\n";

#prepare PATH environmental variable
my $path = File::Basename::dirname ($0);
$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in

$outfile ||= $queryfile;
$genetype ||= 'refgene';
$genetype =~ m/^refgene|knowngene|ensgene|gencodegene$/i or $genetype =~ m/wgEncodeGencode\w+$/ or pod2usage ("Error in argument: the --genetype can be 'refgene', 'knowngene', 'ensgene' or 'gencodegene' only");

if ($genetype eq 'gencodegene') {
	if ($buildver eq 'hg18') {
		$genetype = 'wgEncodeGencodeManualV3';	#this is very outdated now, but kept here for backward compatibility
		print STDERR "WARNING: the --genetype argument is set to wgEncodeGencodeManualV3 by default\n";
	} elsif ($buildver eq 'hg19') {
		$genetype = 'wgEncodeGencodeManualV4';	#this is very outdated now, but kept here for backward compatibility
		print STDERR "WARNING: the --genetype argument is set to wgEncodeGencodeManualV4 by default\n";
	}
}

if (not defined $buildver) {
	$buildver = 'hg18';
	print STDERR "NOTICE: the --buildver argument is set as 'hg18' by default\n";
}

if (defined $maf_threshold) {
	pod2usage ("Error in argument: the --maf_threshold is removed due to user complaints. Please use --aaf_threshold instead");
}


if (defined $aaf_threshold) {
	$aaf_threshold >= 0 and $aaf_threshold <= 1 or pod2usage ("Error: the --aaf_threshold argument must be between 0 and 1 inclusive");
}

not defined $checkfile and $checkfile = 1;

if (not $protocol) {	#kept here for backward compatibility
	$operation and pod2usage ("Error in argument: you must specify --protocol if you specify --operation");
	if ($buildver eq 'hg18') {
		$protocol = 'nonsyn_splicing,1000g2010jul_ceu,1000g2010jul_jptchb,1000g2010jul_yri,snp132,esp5400_ea,esp5400_aa,recessive';
		$operation = 'g,f,f,f,f,f,f,m';
		print STDERR "NOTICE: the --protocol argument is set as 'nonsyn_splicing,1000g2010jul_ceu,1000g2010jul_jptchb,1000g2010jul_yri,snp132,esp5400_ea,esp5400_aa,recessive' by default\n";
	} elsif ($buildver eq 'hg19') {
		$protocol = 'nonsyn_splicing,1000g2012apr_all,snp135NonFlagged,esp6500_ea,esp6500_aa,recessive';
		$operation = 'g,f,f,f,f,m';
		print STDERR "NOTICE: the --protocol argument is set as 'nonsyn_splicing,1000g2012apr_all,snp135NonFlagged,esp6500_ea,esp6500_aa,recessive' by default\n";
	} elsif ($buildver eq 'hg38') {
		$protocol = 'nonsyn_splicing,1000g2015aug_all,snp142,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,recessive';
		$operation = 'g,f,f,f,f,m';
		print STDERR "NOTICE: the --protocol argument is set as 'nonsyn_splicing,1000g2015aug_all,snp142,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,recessive' by default\n";
	}
}

if ($protocol =~ m/\bgeneric\b/) {
	$genericdbfile or pod2usage ("Error in argument: please specify -genericdbfile argument when 'generic' operation is specified");
}

my @protocol = split (/,/, $protocol);
my @operation = split (/,/, $operation);
my @argument = split (/,/, $argument||'', -1);
my $sc;
my $linecount;

my (%valistep, $skip);

@protocol == @operation or pod2usage ("Error in argument: different number of elements are specified in --protocol and --operation argument");
@argument and @protocol == @argument || pod2usage ("Error in argument: different number of elements are specified in --protocol and --argument argument");

for my $op (@operation) {
	$op =~ m/^g|r|f|m|rr|fr$/ or pod2usage ("Error in argument: the --operation argument must be comma-separated list of 'g', 'r', 'f', 'rr', 'fr' or 'm'");
}

$checkfile and checkFileExistence ();

copy ($queryfile, "$outfile.step0.varlist");
for my $i (0 .. @protocol-1) {
	print STDERR "-----------------------------------------------------------------\n";
	print STDERR "NOTICE: Processing operation=$operation[$i] protocol=$protocol[$i]\n";
	if ($operation[$i] eq 'g') {
		geneOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 0, $argument[$i]||undef);
	} elsif ($operation[$i] eq 'r') {
		regionOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 0, $argument[$i]||undef);
	} elsif ($operation[$i] eq 'rr') {
		regionOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 1, $argument[$i]||undef);
	} elsif ($operation[$i] eq 'f') {
		filterOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 0, $argument[$i]||undef);
	} elsif ($operation[$i] eq 'fr') {
		filterOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 1, $argument[$i]||undef);
	} elsif ($operation[$i] eq 'm') {
		modelOperation ($i+1, "$outfile.step$i.varlist", $protocol[$i], 0, $argument[$i]||undef);
	}
}



sub geneOperation {
	my ($step, $infile, $operation, $reverse, $arg) = @_;
	
	if ($operation eq 'nonsyn_splicing' or $operation eq 'nonsyn') {
		$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $genetype -outfile $outfile.step$step $infile $dbloc";
		defined $arg and $sc .= " $arg";
		print STDERR "\nNOTICE: Running step $step with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
		columnGrep ("$outfile.step$step.exonic_variant_function", "$outfile.step$step.varlist", '^(?:synonymous SNV|nonframeshift deletion|nonframeshift insertion|nonframeshift substitution)$', 2, "\t", 4, 'reverse');
		
		if ($operation eq 'nonsyn_splicing') {
			columnGrep ("$outfile.step1.variant_function", "$outfile.step$step.varlist.temp", '\bsplicing\b', 1, "\t", 1);
			columnGrep ("$outfile.step$step.varlist.temp", "$outfile.step$step.varlist", '\bexonic\b', 1, "\t", 3, 'reverse', 'append');
		}

		open (FH, "$outfile.step$step.varlist") or die "Error: cannot read from $outfile.step$step.varlist: $!\n";
		open (FHOUT, ">$outfile.step$step.varlist.temp") or die "Error: cannot write to $outfile.step$step.varlist.temp: $!\n";
		my %found;
		my $linecount=0;
		while (<FH>) {
			$found{$_} and next;
			print FHOUT $_;
			$found{$_}++;
			$linecount++;
		}
		close (FH);
		close (FHOUT);
		move ("$outfile.step$step.varlist.temp", "$outfile.step$step.varlist");
		
		#$remove and unlink ("$outfile.step$step.varlist");
		
		$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
		print STDERR "NOTICE: After step $step, $linecount variants are left in analysis.\n";
	} else {
		die "Error: the $operation command for gene-based annotation is currently not supported\n";
	}
}

sub columnGrep {
	my ($infile, $outfile, $pattern, $column, $separator, $outcolumn, $reverse, $append) = @_;
	defined $column or die "Error: column to the columnGrep subroutine must be defined\n";
	defined $separator or $separator = '\s+';
	open (FHIN, $infile) or die "Error: cannot read from input file $infile: $!\n";
	if ($append) {
		open (FHOUT, ">>$outfile") or die "Error: cannot append to output file $outfile: $!\n";
	} else {
		open (FHOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
	}
	while (<FHIN>) {
		my @field = split (/$separator/, $_);
		if ($reverse) {
			$field[$column-1] =~ m/$pattern/ or print FHOUT join ($separator, @field[($outcolumn-1) .. $#field]);
		} else {
			$field[$column-1] =~ m/$pattern/ and print FHOUT join ($separator, @field[($outcolumn-1) .. $#field]);
		}
	}
	close (FHIN);
	close (FHOUT);
}

sub regionOperation {
	my ($step, $infile, $dbtype, $reverse, $arg) = @_;
	$sc = "annotate_variation.pl -regionanno -dbtype $dbtype -buildver $buildver -outfile $outfile.step$step $infile $dbloc";
	defined $arg and $sc .= " $arg";
	print STDERR "\nNOTICE: Running step $step with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	if ($reverse) {
		system ("cut -f 3- $outfile.step$step.${buildver}_$dbtype > $outfile.step$step.temp");
		system ("fgrep -v -f $outfile.step$step.temp $infile >  $outfile.step$step.varlist");
	} else {
		system ("cut -f 3- $outfile.step$step.${buildver}_$dbtype > $outfile.step$step.varlist");
	}
	
	#$remove and unlink ("$outfile.step$step.varlist", "$outfile.step$step.${buildver}_$dbtype");
	$linecount = qx/cat $outfile.step$step.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
	print STDERR "NOTICE: After step $step, $linecount variants are left in analysis.\n";
}

sub filterOperation {
	my ($step, $infile, $dbtype, $reverse, $arg) = @_;
	$sc = "annotate_variation.pl -filter -dbtype $dbtype -buildver $buildver -outfile $outfile.step$step $infile $dbloc";
	defined $arg and $sc .= " $arg";
	
	if ($dbtype eq 'generic') {
		$sc .= " -genericdbfile $genericdbfile";
	}
	
	if ($reverse) {
		$sc .= " -reverse";
	}
	
	if ($dbtype eq 'ljb_sift') {
		my $score_threshold = $ljb_sift_threshold || 0.95;
		$sc .= " -score_threshold $score_threshold -reverse";
	}
	
	if ($dbtype eq 'ljb_pp2') {
		my $score_threshold = $ljb_pp2_threshold || 0.85;
		$sc .= " -score_threshold $score_threshold -reverse";
	}
	
	if ($dbtype =~ m/ljb2\d*_sift/) {
		my $score_threshold = $ljb2_sift_threshold || 0.05;
		$sc .= " -score_threshold $score_threshold";
	}
	
	if ($dbtype =~ m/ljb2\d*_pp2hdiv/) {
		my $score_threshold = $ljb2_pp2hdiv_threshold || 0.957;
		$sc .= " -score_threshold $score_threshold -reverse";
	}
	
	if ($dbtype =~ m/ljb2\d*_pp2hvar/) {
		my $score_threshold = $ljb2_pp2hvar_threshold || 0.909;
		$sc .= " -score_threshold $score_threshold -reverse";
	}
	
	if (defined $aaf_threshold) {
		if ($dbtype =~ m/^1000g/ or $dbtype =~ m/^esp\d+/ or $dbtype =~ m/^cg\d+/ or $dbtype =~ m/^exac\d+/) {
			$sc .= " -score_threshold $aaf_threshold";
		}
	}
	
	
	print STDERR "\nNOTICE: Running step $step with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	
	my $dbtype1 = $dbtype;
	if ($dbtype =~ m/^1000g_(\w+)/) {
		$dbtype1 = uc ($1) . ".sites.2009_04";
	} elsif ($dbtype =~ m/^1000g2010_(\w+)/) {
		$dbtype1 = uc ($1) . ".sites.2010_03";
	} elsif ($dbtype =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
		my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
		$dbtype1 = uc ($3) . ".sites.$1" . '_' . $monthhash{$2};
	}
		
	copy ("$outfile.step$step.${buildver}_${dbtype1}_filtered", "$outfile.step$step.varlist");		#use dbtype1, not dbtype!!!

	#$remove and unlink ("$outfile.step$step.varlist", "$outfile.step$step.${dbtype}_filtered", "$outfile.step$step.${dbtype}_dropped");
	$linecount = qx/cat $outfile.step$step.varlist | wc -l/; chomp $linecount;
	$linecount or warn "WARNING: No variants were left in analysis after this step. Program exits.\n" and exit;
	print STDERR "NOTICE: After step $step, $linecount variants are left in analysis.\n";
}

sub modelOperation {
	my ($step, $infile, $dbtype, $reverse, $arg) = @_;
	$sc = "fgrep -f $infile $outfile.step1.exonic_variant_function | fgrep -v -w UNKNOWN | cut -f 2- > $outfile.step$step.varlist;";		#function, gene name, plus original input
	$sc .= "cut -f 3- $outfile.step$step.varlist > $outfile.step$step.temp;";			#list of all avinput
	$sc .= "fgrep -v -f $outfile.step$step.temp $infile > $outfile.step$step.temp1;";		#list of splicing variants
	$sc .= "fgrep -f $outfile.step$step.temp1 $outfile.step1.variant_function | fgrep splicing >> $outfile.step$step.varlist;";	#adding splicing variants to nonsyn variants
	print STDERR "\nNOTICE: Running step 8 with system command <$sc>\n";
	system ($sc);			#this command may generate error, because the $outfile.step8.temp1 file may be empty

	$remove and unlink ("$outfile.step$step.temp", "$outfile.step$step.temp1", "$outfile.step1.exonic_variant_function");


	my (%found, %varpos);		#count of gene, variant information of the variant
	open (VAR, "$outfile.step$step.varlist") or die "Error: cannot read from varlist file $outfile.step$step.varlist: $!\n";
	while (<VAR>) {
		my @field = split (/\t/, $_);
		$field[1] =~ s/,$//;
		
		#$field[1] =~ s/\([^\(\)]+\)//g;		#handle situations such as splicing        EMG1(NM_006331:exon1:c.125+1T>GC,NM_006331:exon2:c.126-1T>GC)   
		
		$field[1] =~ m/^(\w+)/ or die "Error: invalid record in input file $outfile.step$step.varlist (gene name expected at second column): <$_>\n";
		my $gene = $1;
		$found{$gene}++;
		$varpos{$gene} .= "\t$field[1]";
		if (m/\bhom\b/) {
			$found{$gene}++;
		}
	}
	
	my $count_candidate_gene = 0;
	open (OUT, ">$outfile.step$step.genelist") or die "Error: cannot write to output file $outfile.step$step.genelist: $!\n";
	print OUT "Gene\tNumber_of_deleterious_alleles\tMutations\n";
	for my $key (keys %found) {
		if ($dbtype eq 'recessive') {
			if ($found{$key} >= 2) {
				print OUT "$key\t$found{$key}$varpos{$key}\n";
				$count_candidate_gene++;
			}
		} elsif ($dbtype eq 'dominant') {
			if ($found{$key} >= 1) {
				print OUT "$key\t$found{$key}$varpos{$key}\n";
				$count_candidate_gene++;
			}
		} else {
			die "Error: the model operation $dbtype specified in -operation argument is not supported\n";
		}
	}
	print STDERR "\nNOTICE: a list of $count_candidate_gene potentially important genes and the number of deleterious alleles in them are written to $outfile.step$step.genelist\n";
}

sub checkFileExistence {
	my @file;
	my %dbtype1 = ('gene'=>'refGene', 'refgene'=>'refGene', 'knowngene'=>'knownGene', 'ensgene'=>'ensGene', 'band'=>'cytoBand', 'cytoband'=>'cytoBand', 'tfbs'=>'tfbsConsSites', 'mirna'=>'wgRna',
			'mirnatarget'=>'targetScanS', 'segdup'=>'genomicSuperDups', 'omimgene'=>'omimGene', 'gwascatalog'=>'gwasCatalog', 
			'1000g_ceu'=>'CEU.sites.2009_04', '1000g_yri'=>'YRI.sites.2009_04', '1000g_jptchb'=>'JPTCHB.sites.2009_04', 
			'1000g2010_ceu'=>'CEU.sites.2010_03', '1000g2010_yri'=>'YRI.sites.2010_03', '1000g2010_jptchb'=>'JPTCHB.sites.2010_03',
			'1000g2010jul_ceu'=>'CEU.sites.2010_07', '1000g2010jul_yri'=>'YRI.sites.2010_07', '1000g2010jul_jptchb'=>'JPTCHB.sites.2010_07',
			'1000g2010nov_all'=>'ALL.sites.2010_11', '1000g2011may_all'=>'ALL.sites.2011_05'
			);		#for backward compatibility
	for my $i (0 .. @protocol-1) {
		my $dbtype1;
		if ($operation[$i] eq 'g') {
			$dbtype1 = $dbtype1{$genetype} || $genetype;
		} elsif ($operation[$i] eq 'm') {
			next;
		} else {
			$dbtype1 = $dbtype1{$protocol[$i]} || $protocol[$i];
		}
		
		if ($protocol[$i] =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$dbtype1 = uc ($3) . '.sites.' . $1 . '_' . $monthhash{$2};
		}
		my $file;
		if ($dbtype1 ne 'generic') {
			$file = $buildver . "_" . $dbtype1 . ".txt";
			push @file, $file;
		}
	}

	for my $i (0 .. @file-1) {
		my $dbfile = File::Spec->catfile ($dbloc, $file[$i]);
		-f $dbfile or die "Error: the required database file $dbfile does not exist. Please download it via -downdb argument by annotate_variation.pl.\n";
	}
}


=head1 SYNOPSIS

 variants_reduction.pl [arguments] <query-file> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --protocol <string>		comma-delimited string specifying annotation protocol
            --operation <string>	comma-delimited string specifying type of operation
            --outfile <string>		output file name prefix
            --buildver <string>		genome build version (default: hg18)
            --remove			remove all temporary files
            --genetype <string>		gene definition (default: refgene)
            --aaf_threshold <float>	alternative allele frequency threshold for filtering
            --(no)checkfile		check if database file exists (default: ON)
            --genericdbfile <file>	specify generic db file
            --ljb_sift_threshold <float>	specify the threshold for ljb_sift (default: 0.95)
            --ljb_pp2_threshold <float>		specify the threshold for ljb_pp2 (default: 0.85)
            --ljb2_sift_threshold <float>	specify the threshold for ljb2_sift (default: 0.05)
            --ljb_pp2hvar_threshold <float>	specify the threshold for ljb2_pp2hvar (default: 0.909)
            --ljb_pp2hdiv_threshold <float>	specify the threshold for ljb2_pp2hvar (default: 0.957)
            --argument <string>		commad-delimited strings as optional argument for each operation
            

 Function: automatically run a pipeline on a list of variants (potentially 
 whole-genome SNPs from a patient with Mendelian disease) and identify a small 
 subset that are most likely causal for Mendelian diseases
 
 Example: variants_reduction.pl sample.avinput humandb/ -protocol nonsyn_splicing,genomicSuperDups,phastConsElements46way,1000g2012apr_all,esp5400_ea,esp5400_aa,snp135NonFlagged,dominant -operation g,rr,r,f,f,f,f,m -out reduce -buildver hg19
                  
 Version: $Date: 2015-12-14 13:51:19 -0800 (Mon, 14 Dec 2015) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--protocol>

comma-delimited string specifying annotation protocol. These strings typically 
represent database names in ANNOVAR.

=item B<--operation>

comma-delimited string specifying type of operation. These strings can be g 
(gene), r (region), rr (region reverse), f (filter) or fr (filter reverse). The 
reverse suffix is needed for scores where higher ones are functionally important 
ones, unlike other scores (such as allele frequency) where smaller scores are 
functionally important ones. For ljb_sift/ljb2_sift and ljb_pp2/ljb2_pp2, the 
reverse mechanism is set automatically for users to make life easier.

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--genetype>

specify the gene definition, such as refgene (default), ucsc known gene, ensembl 
gene and gencode gene.

=item B<--aaf_threshold>

specify the alternative allele frequency threshold. This argument works 
for 1000 Genomes Project, ESP database and CG (complete genomics) database.

=item B<--(no)checkfile>

the program will check if all required database files exist before execution of annotation

=item B<--genericdbfile>

specify the genericdb file used in -dbtype generic

=item B<--ljb_sift_threshold>

specify the LJB_SIFT threshold for filter operation (default: 0.95). NOTE THAT 
IN LJB_SIFT, THE SIFT SCORE IS REPORTED AS 1-SIFT.

=item B<--ljb_pp2_threshold>

specify the LJB_PP2 threshold for filter operation (default: 0.85)

=item B<--ljb2_sift_threshold>

specify the LJB2_SIFT threshold for filter operation (default: 0.05). NOTE THAT 
IN LJB2_SIFT, SIFT SCORE IS REPORTED AS THE ORIGINAL SCORE.

=item B<--ljb2_pp2hdiv_threshold>

specify the LJB_PP2DIV threshold for filter operation (default: 0.957)

=item B<--ljb_pp2hvar_threshold>

specify the LJB_PP2_HVAR threshold for filter operation (default: 0.909)

=item B<--argument>

comma-delimited argument list to be supplied to each of the operations. This 
allows users to fine-tune the reduction procedure.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, such as those generated from next-generation sequencing 
experiments. 

The variants_reduction.pl program within the ANNOVAR package 
provides users with the ability to perform stepwise variants reduction on a 
large set of input variants, to help trim down to a list of functionally 
important variants. Through the --protocol argument, users can specify a list of 
procedures to filter down variants, step by step.

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact $Author: Kai Wang <kai@openbioinformatics.org> $.

=cut