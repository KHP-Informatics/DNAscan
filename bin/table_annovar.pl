#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use File::Spec;

our $REVISION = '$Revision: 11f4bb33cd1289ba7c43dadc20b3e9982b5a9a00 $';
our $DATE =	'$Date: 2016-02-01 00:11:18 -0800 (Mon,  1 Feb 2016) $';  
our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org> $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $remove, $checkfile, $protocol, $operation, $otherinfo, $onetranscript, $nastring, $genericdbfile, $gff3dbfile, $bedfile, $vcfdbfile, $csvout, $argument, $tempdir, $vcfinput, $dot2underline, $thread, $maxgenethread);
our $orig_command = $0;

for my $i (0 .. @ARGV-1) {
	my $temp = $ARGV[$i];
	if ($temp =~ m/\s/) {
		my @temp = split (/,/, $temp, -1);
		$temp = join (',', map {qq#'$_'#} @temp);	#handle situations like -argument '--hgvs --exonicsplicing',,,,,,,,,,,,,,,, in command line
	}
	$orig_command .= " $temp";
	#print STDERR "orig_command = $orig_command\n";
}

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'remove'=>\$remove, 'checkfile!'=>\$checkfile, 
	'protocol=s'=>\$protocol, 'operation=s'=>\$operation, 'otherinfo'=>\$otherinfo, 'onetranscript'=>\$onetranscript, 'nastring=s'=>\$nastring,
	'genericdbfile=s'=>\$genericdbfile, 'gff3dbfile=s'=>\$gff3dbfile, 'bedfile=s'=>\$bedfile, 'vcfdbfile=s'=>\$vcfdbfile,
	'csvout'=>\$csvout, 'argument=s'=>\$argument, 'tempdir=s'=>\$tempdir, 'vcfinput'=>\$vcfinput, 'dot2underline'=>\$dot2underline,
	'thread=i'=>\$thread, 'maxgenethread=i'=>\$maxgenethread) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($queryfile, $dbloc) = @ARGV;


my @unlink;					#a list of files to be deleted
my @header;					#annotation output header
my %varanno;					#varstring as key of 1st hash, anno_db as key of 2nd hash, anno as value of 2nd hash
my (@protocol, @operation, @dbtype1, @argument);		#@dbtype1 is the translated file name for @protocol
my (@genericdbfile, @gff3dbfile, @vcfdbfile, @bedfile);
my (@genericdbfile1, @gff3dbfile1, @vcfdbfile1, @bedfile1);
my ($tempfile);			#prefix for temporary files (by default, temporary files are written to the -outfile directory

my %annotation_headers = (
	"ljb_all" => [ qw/LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++/ ],
	"ljb2_all" => [ qw/LJB2_SIFT LJB2_PolyPhen2_HDIV LJB2_PP2_HDIV_Pred LJB2_PolyPhen2_HVAR LJB2_PolyPhen2_HVAR_Pred LJB2_LRT LJB2_LRT_Pred LJB2_MutationTaster LJB2_MutationTaster_Pred LJB_MutationAssessor LJB_MutationAssessor_Pred LJB2_FATHMM LJB2_GERP++ LJB2_PhyloP LJB2_SiPhy/ ],
	"ljb23_all" => [ qw/LJB23_SIFT_score LJB23_SIFT_score_converted LJB23_SIFT_pred LJB23_Polyphen2_HDIV_score LJB23_Polyphen2_HDIV_pred LJB23_Polyphen2_HVAR_score LJB23_Polyphen2_HVAR_pred LJB23_LRT_score LJB23_LRT_score_converted LJB23_LRT_pred LJB23_MutationTaster_score LJB23_MutationTaster_score_converted LJB23_MutationTaster_pred LJB23_MutationAssessor_score LJB23_MutationAssessor_score_converted LJB23_MutationAssessor_pred LJB23_FATHMM_score LJB23_FATHMM_score_converted LJB23_FATHMM_pred LJB23_RadialSVM_score LJB23_RadialSVM_score_converted LJB23_RadialSVM_pred LJB23_LR_score LJB23_LR_pred LJB23_GERP++ LJB23_PhyloP LJB23_SiPhy/ ],
	"popfreq_all" => [ qw/PopFreqMax 1000G2012APR_ALL 1000G2012APR_AFR 1000G2012APR_AMR 1000G2012APR_ASN 1000G2012APR_EUR ESP6500si_ALL ESP6500si_AA ESP6500si_EA CG46/ ]
	);		#for backward compatibility (older annotation databases do not have the comment line)

processArgument ();
@dbtype1 = proxyDBType(@protocol);
$checkfile and checkFileExistence (@dbtype1);

if ($vcfinput) {
	my $sc;
	$csvout and pod2usage ("Error in argument: -csvout is not compatible with -vcfinput");
	$sc = "convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 $queryfile > $tempfile.avinput";
	print STDERR "\nNOTICE: Running with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	$sc = $orig_command;
	$sc =~ s/\s$queryfile\s/ $tempfile.avinput /;		#change command line argument
	$sc =~ s/\s-?-vcfi\w*//;				#delete -vcfinput argument
	$sc .= ' -otherinfo';					#add -otherinfo so that VCF information is included in output file
	$sc =~ m/\s\-na/ or $sc .= ' -nastring .';		#force the nastring to be dot in output VCF file
	if ($sc =~ s/\s-?-out\w*\s+\S+/ -outfile $tempfile/) {	#change the -outfile argument to the temporary location
		1;
	} else {
		$sc .= " -outfile $tempfile";
	}
	print STDERR "\nNOTICE: Running with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	#push @unlink, "$tempfile.${buildver}_multianno.txt";	#by default, we no longer delete the tab-delimited file
	backConvertVCF ("$tempfile.${buildver}_multianno.txt", $queryfile);
} else {
	for my $i (0 .. @protocol-1) {
		print STDERR "-----------------------------------------------------------------\n";
		print STDERR "NOTICE: Processing operation=$operation[$i] protocol=$protocol[$i]\n";
		if ($operation[$i] eq 'g') {
			geneOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
		} elsif ($operation[$i] eq 'r') {
			regionOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
		} elsif ($operation[$i] eq 'f') {
			filterOperation ($protocol[$i], $dbtype1[$i], $argument[$i]||undef);
		}
	}
	
	printOrigOutput ();
}

# Remove temporary files when -remove argument is specified
if ($remove) {
	unlink(@unlink);
	if ($tempdir) {
		$tempfile =~ s/temp$//;
		rmdir $tempfile;
	}
}

# Given the input VCF file and the tab-delimited annotation file, generate a new VCF that contains the annotation in the INFO column
sub backConvertVCF {
	my ($multiannofile, $vcfin) = @_;
	open (MANNO, $multiannofile) or die "Error: cannot read from multianno file $multiannofile: $!\n";
	print STDERR "NOTICE: Reading from $multiannofile\n";
	$_ = <MANNO>;
	s/[\r\n]+$//;
	my @name = split (/\t/, $_);
	$name[$#name] eq 'Otherinfo' or die "Error: the last column in header row should be 'Otherinfo'\n";
	
	my @nextrecord;
	my (@field, @prefield);
	my ($anno_string, $pre_anno_string);
	
	open (OUT, ">$outfile.${buildver}_multianno.vcf") or die "Error: cannot write to output file $outfile.${buildver}_multianno.vcf\n";
	print STDERR "-----------------------------------------------------------------\n";
	print STDERR "NOTICE: VCF output is written to $outfile.${buildver}_multianno.vcf\n";

	#open (IN, $vcfin) or die "Error: cannot read from input VCF file $vcfin: $!\n";

	if ($vcfin =~ m/\.gz$/) {		#handle vcf.gz file
		open (IN, "gunzip -c $vcfin |") or die "Error: cannot read from input VCF file $vcfin: $!\n";
	} else {
		open (IN, $vcfin) or die "Error: cannot read from input VCF file $vcfin: $!\n";
	}

	while (<IN>) {
		s/[\r\n]+$//;
		if (m/^##/) {
			print OUT $_, "\n";
		} elsif (m/^#CHROM/) {
			print OUT qq{##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,Description="Flag the start of ANNOVAR annotation for one alternative allele">\n};
			for my $i (0 .. @name-1) {
				if ($name[$i] =~ m/^(Chr|Start|End|Ref|Alt)$/) {
					1;					#these are annovar input and are not in INFO field when printed out
				} elsif ($name[$i] =~ m/^(1000g\d+|esp\d+|cg\d+|popfreq|nci\d+)/) {
					print OUT qq{##INFO=<ID=$name[$i],Number=1,Type=Float,Description="$name[$i] annotation provided by ANNOVAR">\n};
				} elsif ($name[$i] eq 'Otherinfo') {		#this is the field that stores VCF information, and will not be in the annotation in INFO field
					1;
				} else {
					print OUT qq{##INFO=<ID=$name[$i],Number=.,Type=String,Description="$name[$i] annotation provided by ANNOVAR">\n};
				}
			}
			print OUT qq{##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description="Flag the end of ANNOVAR annotation for one alternative allele">\n};
			print OUT $_, "\n";
			last;
		} elsif (m/^#/) {
			print OUT $_, "\n";
		} else {
			last;		#could be a VCF file without a valid head line
		}
	}
	close (IN);
	
	while (<MANNO>) {
		s/[\r\n]+$//;
		my $anno_string;
		@field = split (/\t/, $_);
		
		$DATE =~ m/Date: (\d+\-\d+\-\d+)/;
		$anno_string = ";ANNOVAR_DATE=$1";
		for my $i (5 .. @name-2) {		#these are all the annotation columns (starting from 6th column to the last one, which is Otherinfo which will not be written)
			$field[$i] =~ s/\s/_/g;		#convert 'nonsynonymous SNV' to 'nonsynonymous' as space is not allowed in VCF INFO field, yet _ is easier to view for users
			$field[$i] =~ s/;/\\x3b/g;	#; is not allowed in VCF INFO field
			$field[$i] =~ s/=/\\x3d/g;	#= is not allowed in VCF INFO field
			$anno_string .= ";$name[$i]=$field[$i]";	#the 8th column in Otherinfo is INFO column (plus the freq, quality, dp column)
		}
		$anno_string .= ";ALLELE_END";
		
		if (@prefield and join ("\t", @field[@name+3 .. $#field]) eq join ("\t", @prefield[@name+3 .. $#field])) {		#same locus, multiple alternative allele
			$pre_anno_string .= $anno_string;
			
		} else {
			if (@prefield) {
				$prefield[@name+9] .= $pre_anno_string;		#update INFO field
				print OUT join ("\t", @prefield[@name+2 .. $#prefield]), "\n";
			}
				
			@prefield = @field;
			$pre_anno_string = $anno_string;
		}
	}
	$prefield[@name+9] .= $pre_anno_string;
	print OUT join ("\t", @prefield[@name+2 .. $#prefield]), "\n";
}

# Print out the annotation using the same order as the input file (previous version of table_annovar has random order which some users complained about; however, now the output will contain duplicated lines if user have duplicated variants in input).
sub printOrigOutput {
	my $final_out;
	if ($csvout) {
		$final_out="$outfile.${buildver}_multianno.csv";
	} else {
		$final_out="$outfile.${buildver}_multianno.txt";
	}
	open OUT,">",$final_out or die "Cannot write to $final_out: $!\n";
	print STDERR "-----------------------------------------------------------------\n";
	print STDERR "NOTICE: Multianno output file is written to $final_out\n";
	
	my @expanded_header;

	for my $item (@header) {
		if ( exists $annotation_headers{$item} ) {
			push @expanded_header, @{ $annotation_headers{$item} };
		} else {
			push @expanded_header, $item;
		}
	}
	
	if ($csvout) {
		print OUT join(",", qw/Chr Start End Ref Alt/, @expanded_header), $otherinfo?",Otherinfo":"", "\n";
	} else {
		print OUT join("\t", qw/Chr Start End Ref Alt/, @expanded_header), $otherinfo?"\tOtherinfo":"", "\n";
	}
	
	open (FH, $queryfile) or die "Error: cannot read from inputfile $queryfile: $!\n";
	while (<FH>) {
		s/[\r\n]+$//;
		m/^(\S+\s+\S+\s+\S+\s+\S+\s+\S+)\s*(.*)/ or next;
		my ($varstring, $info) = ($1, $2);
		my @varstring = split (/\s+/, $varstring);
		$varstring =~ s/\s+/\t/g;
		
		my @oneline;
		for my $i (0 .. @header-1) {
			my $item = $header[$i];
			my $expanded_field;

			if ($dot2underline) {
				$item =~ s/Func_/Func./;
				$item =~ s/Gene_/Gene./;
				$item =~ s/GeneDetail_/GeneDetail./;
				$item =~ s/ExonicFunc_/ExonicFunc./;
				$item =~ s/AAChange_/AAChange./;
			}
			
			if ( exists $annotation_headers{$item} ) {
				$expanded_field = scalar @{ $annotation_headers{$item} };
			}
			
			if ($csvout) {
				if ($expanded_field) {
					push @oneline,($varanno{$varstring}{$item} || join (",", ($nastring) x $expanded_field));	#value is already comma-delimited
				} else {
					push @oneline,(defined $varanno{$varstring}{$item}) ? qq/"$varanno{$varstring}{$item}"/ : $nastring;	#value of $varanno{$varstring}{$item} could be zero
				}
			} else {
				if ($expanded_field) {
					if ($varanno{$varstring}{$item}) {
						for my $nextscore (split (/,/, $varanno{$varstring}{$item})) {
							$nextscore =~ s/\\x2c/,/g;	#now it is safe to change \x2c to , (after split operation above)
							$nextscore =~ s/\\x23/#/g;	#now it is safe to change \x23 to # (after split operation above)
							push @oneline, $nextscore;
						}
					} else {
						for (1 .. $expanded_field) {
							push @oneline, $nastring;
						}
					}
				} else {
					push @oneline, (defined $varanno{$varstring}{$item})? $varanno{$varstring}{$item} : $nastring;	#value of $varanno{$varstring}{$item} could be zero
				}
			}
		}
		if ($csvout) {
			print OUT join (",", @varstring, @oneline), $otherinfo?qq/,"$info"/:"", "\n";
		} else {
			print OUT join ("\t", @varstring, @oneline), $otherinfo?"\t$info":"", "\n";
		}
	}
}

# Process the arguments supplied by users
sub processArgument {
	$outfile ||= $queryfile;

	if ($tempdir) {
		$tempdir =~ s/[\\\/]+$//;		#remove trailing slash
		my $maxLenth=8;
		my @a = (0..9,'a'..'z','A'..'Z','-','_');
		my $password;
		while (1) {
			$password = join '', map { $a[int rand @a] } 0..($maxLenth-1);
			my $tempsubdir = File::Spec->catfile ($tempdir, $password);
			-d $tempsubdir and next;
			mkdir $tempsubdir or die "Error: -tempdir ($tempdir) is not writtable\n";
			print STDERR "NOTICE: temporary files will be written to $tempsubdir\n";
			$tempfile = File::Spec->catfile ($tempsubdir, 'temp');
			last;
		}
	} else {
		$tempfile = $outfile;
	}
        	
	if (not defined $buildver) {
		$buildver = 'hg18';
		print STDERR "NOTICE: the --buildver argument is set as 'hg18' by default\n";
	}
	
	not defined $checkfile and $checkfile = 1;
	
	
	if ($vcfinput) {
		defined $nastring and $nastring ne '.' and pod2usage ("Error in argument: -nastring must be '.' when '-vcfinput' is specified");
		$nastring = '.';		#force "." to denote missing data. "." can be recognized as NA in other tools such as R
	} else {
		defined $nastring or $nastring = '';
	}
	
	if (not $protocol) {
		$operation and pod2usage ("Error in argument: you must specify --protocol if you specify --operation");
		if ($buildver eq 'hg18') {
			$protocol = 'gene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb_all';
			$operation = 'g,r,r,f,f,f,f';
			print STDERR "NOTICE: the --protocol argument is set as 'gene,phastConsElements44way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all' by default\n";
		} elsif ($buildver eq 'hg19') {
			$protocol = 'gene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb_all';
			$operation = 'g,r,r,f,f,f,f';
			print STDERR "NOTICE: the --protocol argument is set as 'gene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp137,ljb2_all' by default\n";
		} else {
			pod2usage ("Error in argument: please specify --protocol argument for the --buildver $buildver");
		}
	}
	
	if ($protocol =~ m/\bgeneric\b/) {
		$genericdbfile or pod2usage ("Error in argument: please specify -genericdbfile argument when 'generic' operation is specified");
	}
	
	#additional preparation work
	@protocol = split (/,/, $protocol);
	@operation = split (/,/, $operation);
	@argument = split (/,/, $argument||'', -1);
	@protocol == @operation or pod2usage ("Error in argument: different number of elements are specified in --protocol and --operation argument");
	@argument and @protocol == @argument || pod2usage ("Error in argument: different number of elements are specified in --protocol and --argument argument");
	for my $op (@operation) {
		$op =~ m/^g|r|f$/ or pod2usage ("Error in argument: the --operation argument must be comma-separated list of 'g', 'r', 'f'");
	}
		
	my %uniq_protocol;
	for (@protocol) {
		$uniq_protocol{$_}++;
	}
	if ($uniq_protocol{generic}) {
		$genericdbfile or pod2usage ("Error in argument: please specify -genericdbfile argument when 'generic' protocol is specified");
		@genericdbfile = split (/,/, $genericdbfile);
		@genericdbfile == $uniq_protocol{generic} or pod2usage ("Error in argument: you specified $uniq_protocol{generic} 'generic' in 'protocol' argument, but only ${\(scalar @genericdbfile)} filenames in 'genericdbfile' argument");
	}
	if ($uniq_protocol{gff3}) {
		$gff3dbfile or pod2usage ("Error in argument: please specify -gff3dbfile argument when 'gff3' protocol is specified");
		@gff3dbfile = split (/,/, $gff3dbfile);
		@gff3dbfile == $uniq_protocol{gff3} or pod2usage ("Error in argument: you specified $uniq_protocol{gff3} 'gff3' in 'protocol' argument, but only ${\(scalar @gff3dbfile)} filenames in 'gff3dbfile' argument");
	}
	if ($uniq_protocol{vcf}) {
		$vcfdbfile or pod2usage ("Error in argument: please specify -vcfdbfile argument when 'vcf' protocol is specified");
		@vcfdbfile = split (/,/, $vcfdbfile);
		@vcfdbfile == $uniq_protocol{vcf} or pod2usage ("Error in argument: you specified $uniq_protocol{vcf} 'vcf' in 'protocol' argument, but only ${\(scalar @vcfdbfile)} filenames in 'vcfdbfile' argument");
	}
	if ($uniq_protocol{bed}) {
		$bedfile or pod2usage ("Error in argument: please specify -bedfile argument when 'bed' protocol is specified");
		@bedfile = split (/,/, $bedfile);
		@bedfile == $uniq_protocol{bed} or pod2usage ("Error in argument: you specified $uniq_protocol{bed} 'bed' in 'protocol' argument, but only ${\(scalar @bedfile)} filenames in 'bedfile' argument");
	}

	if (defined $thread) {
		$thread > 0 or pod2usage ("Error: the --thread argument must be a positive integer");
	}

	#prepare PATH environmental variable
	my $path = File::Basename::dirname ($0);
	$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in
}

# Call ANNOVAR for gene annotations
sub geneOperation {
	my ($protocol, $dbtype1, $argument) = @_;
	#generate gene anno
	my $sc;
	$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $protocol -outfile $tempfile.$protocol -exonsort $queryfile $dbloc";
	$argument and $sc .= " $argument";
	
	if ($thread) {
		$sc .= " -thread $thread";
	}
	if ($maxgenethread) {
		$sc .= " -maxgenethread $maxgenethread";
	}
	
	print STDERR "\nNOTICE: Running with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
	
	#read in gene anno
	my $anno_outfile="$tempfile.$protocol.variant_function";
	my $e_anno_outfile="$tempfile.$protocol.exonic_variant_function";
	
	open (FUNCTION, "<",$anno_outfile) or die "Error: cannot read from $anno_outfile: $!\n";	
	open (EFUNCTION,"<",$e_anno_outfile) or die "Error: cannot read from $e_anno_outfile: $!\n";
	
	if ($dot2underline) {
		push @header,"Func_$protocol", "Gene_$protocol", "GeneDetail_$protocol", "ExonicFunc_$protocol", "AAChange_$protocol"; #header
	} else {
		push @header,"Func.$protocol", "Gene.$protocol", "GeneDetail.$protocol", "ExonicFunc.$protocol", "AAChange.$protocol"; #header
	}
	
	while (<FUNCTION>) 
	{
		s/[\r\n]+$//;
		m/^([^\t]+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
		my ($function, $gene, $varstring) = ($1, $2, $3);
		my @spliceanno;
		$varstring =~ s/\s+/\t/g;
		#$varanno{$varstring}{"Func.$protocol"}=$function;	#changed to below to handle -arg '-separate' when users want to see output from -separate argument 20140711
		$varanno{$varstring}{"Func.$protocol"}=exists $varanno{$varstring}{"Func.$protocol"} ? $varanno{$varstring}{"Func.$protocol"} . ";$function" : $function;
		
		while ($gene =~ m/\((.+?)\)/g) {
			$gene =~ s/\((.+?)\)//;
			push @spliceanno, $1;
		}
		$varanno{$varstring}{"Gene.$protocol"} = $gene;
		$varanno{$varstring}{"GeneDetail.$protocol"} = join (";", @spliceanno) || $nastring;
	}
	close FUNCTION;
	while (<EFUNCTION>)
	{
		m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
		my ($efunc, $aachange, $varstring) = ($1, $2, $3);
		
		$aachange =~ s/,$//;	#delete the trailing comma
		
		$varstring =~ s/\s+/\t/g;
		my @aachange = split (/:|,/, $aachange);

		#$varanno{$varstring}{"ExonicFunc.$protocol"}=$efunc;	#changed to below to handle -arg '-separate' when users want to see output from -separate argument 20140711
		$varanno{$varstring}{"ExonicFunc.$protocol"}=exists $varanno{$varstring}{"ExonicFunc.$protocol"} ? $varanno{$varstring}{"ExonicFunc.$protocol"} . ";$efunc" : $efunc;
		if (not $onetranscript) {
			#$varanno{$varstring}{"AAChange.$protocol"}=$aachange;	#changed to below to handle -arg '-separate' when users want to see output from -separate argument 20140711
			$varanno{$varstring}{"AAChange.$protocol"}=exists $varanno{$varstring}{"AAChange.$protocol"} ? $varanno{$varstring}{"AAChange.$protocol"} . ";$aachange" : $aachange;
		} else {
			if (@aachange >= 5) 
			{
			    $varanno{$varstring}{"AAChange.$protocol"}="$aachange[1]:$aachange[3]:$aachange[4]"; #only output aachange in first transcript
			} else {
			    $varanno{$varstring}{"AAChange.$protocol"}=$aachange;		#aachange could be "UNKNOWN"
			}
		}
	}
	close EFUNCTION;
	push @unlink, $anno_outfile, $e_anno_outfile, "$tempfile.$protocol.log";
}

# Call ANNOVAR for region annotations
sub regionOperation {
	my ($dbtype, $dbtype1, $argument) = @_;
	my ($userfile, $header);
	my $sc = "annotate_variation.pl -regionanno -dbtype $dbtype -buildver $buildver -outfile $tempfile $queryfile $dbloc";
	$argument and $sc .= " $argument";
	
	if ($thread) {
		$sc .= " -thread $thread";
	}
	
	if ($dbtype eq 'bed') {
		$userfile = shift @bedfile1;
		$sc .= " -bedfile $userfile";
	} elsif ($dbtype eq 'gff3') {
		$userfile = shift @gff3dbfile1;
		$sc .= " -gff3dbfile $userfile";
	}
	
	print STDERR "\nNOTICE: Running with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";

	if ($dbtype eq 'bed') {
		if (@bedfile-@bedfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@bedfile-@bedfile1);
		}
	} elsif ($dbtype eq 'gff3') {
		if (@gff3dbfile-@gff3dbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@gff3dbfile-@gff3dbfile1);
		}
	} else {
		$header = $dbtype;
	}
	push @header,$header;
	
	open (FH, "$tempfile.${buildver}_$dbtype1") or die "Error: cannot open file\n";
	while (<FH>) 
	{
		m/^([^\t]+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
		my ($db,$anno,$varstring)=($1,$2,$3);
		$varstring =~ s/\s+/\t/g;
		#preprocess $anno
		$varanno{$varstring}{$header}=$anno;
	}
	close (FH);
	
	push @unlink, "$tempfile.${buildver}_$dbtype1", "$tempfile.log";
}

# Call ANNOVAR for filter annotations
sub filterOperation {
	my ($dbtype, $dbtype1, $argument) = @_;
	my ($userfile, $header);
	my $sc = "annotate_variation.pl -filter -dbtype $dbtype -buildver $buildver -outfile $tempfile $queryfile $dbloc";
	$argument and $sc .= " $argument";
	if ($thread) {
		$sc .= " -thread $thread";
	}
	
	if ($dbtype eq 'generic') {
		$userfile = shift @genericdbfile1;
		$sc .= " -genericdbfile $userfile";
	} elsif ($dbtype eq 'vcf') {
		$userfile = shift @vcfdbfile1;
		$sc .= " -vcfdbfile $userfile";
	} else {
		my $dbfile = File::Spec->catfile ($dbloc, $buildver . '_' . $dbtype . ".txt");

		if (open (DB, $dbfile)) {	#the first line may start with # which is a header line
			$_ = <DB>;
			s/[\r\n]+$//;
			if (m/^#/) {
				my @field = split (/\t/, $_);
				splice (@field, 0, 5);			#remove the first five columns (the rest are column names
				$annotation_headers{$dbtype} = [@field];
				print STDERR "NOTICE: Finished reading ", scalar (@field), " column headers for '-dbtype $dbtype'\n";
				$sc .= " -otherinfo";			#20150322: if header information is available, we will automatically add '-otherinfo' argument in the command line to ensure that all columns are included in the output file
			} elsif ($dbtype =~ m/^ljb\d*/ or $dbtype =~ m/^popfreq/ or $dbtype =~ m/^custom/) {		#this is for backward compatibility for databases that do not have column headers (although these db are obselete, for various reasons many labs must stick with old versions)
				$sc .= " -otherinfo";
			}
			close (DB);
		}
	}

	if ($dbtype eq 'avsift') {
		$sc .= " -sift_threshold 0";		#for historical reasons (as the -sift_threshold 0.05 was added automatically for avsift historically)
	}
	
	print STDERR "\nNOTICE: Running system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";

	if ($dbtype eq 'generic') {
		if (@genericdbfile-@genericdbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@genericdbfile-@genericdbfile1);
		}
	} elsif ($dbtype eq 'vcf') {
		if (@vcfdbfile-@vcfdbfile1 == 1) {
			$header = $dbtype;
		} else {
			$header = $dbtype . (@vcfdbfile-@vcfdbfile1);
		}
	} else {
		$header = $dbtype;
	}
	push @header,$header;
	
	open (FH, "$tempfile.${buildver}_${dbtype1}_dropped") or die "Error: cannot open file\n";
	while (<FH>) 
	{
		m/^([^\t]+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+).*/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
		my ($db,$anno,$varstring)=($1,$2,$3);
		#$anno =~ s/\\x2c/,/g;			#when annotating filter database with multiple annotations, the comma is replaced by \x23, when -otherinfo is used (#later I commented out, since we cannot do it here, otherwise extra columns are printed)
		#$anno =~ s/\\x23/#/g;			#when annotating filter database with multiple annotations, the comma is replaced by \x2c, when -otherinfo is used
		$varstring =~ s/\s+/\t/g;
		$varanno{$varstring}{$header}=$anno;
	}
	close (FH);
	
	push @unlink, "$tempfile.${buildver}_${dbtype1}_dropped", "$tempfile.${buildver}_${dbtype1}_filtered", "$tempfile.log";
}

# Generate alias for DB. This is for backward compatibilit only.
sub proxyDBType {
	my @db_names = @_;

	my %dbalias=(					#for backward compatibility (historical reasons)
		'gene'=>'refGene', 
		'refgene'=>'refGene', 
		'knowngene'=>'knownGene', 
		'ensgene'=>'ensGene', 
		'band'=>'cytoBand', 
		'cytoband'=>'cytoBand', 
		'tfbs'=>'tfbsConsSites', 
		'mirna'=>'wgRna',
		'mirnatarget'=>'targetScanS', 
		'segdup'=>'genomicSuperDups', 
		'omimgene'=>'omimGene', 
		'gwascatalog'=>'gwasCatalog',
		'1000g_ceu'=>'CEU.sites.2009_04', 
		'1000g_yri'=>'YRI.sites.2009_04', 
		'1000g_jptchb'=>'JPTCHB.sites.2009_04',
		'1000g2010_ceu'=>'CEU.sites.2010_03', 
		'1000g2010_yri'=>'YRI.sites.2010_03', 
		'1000g2010_jptchb'=>'JPTCHB.sites.2010_03',
		'1000g2010jul_ceu'=>'CEU.sites.2010_07', 
		'1000g2010jul_yri'=>'YRI.sites.2010_07', 
		'1000g2010jul_jptchb'=>'JPTCHB.sites.2010_07',
		'1000g2010nov_all'=>'ALL.sites.2010_11', 
		'1000g2011may_all'=>'ALL.sites.2011_05'
	);
    
	for (@db_names) {
		$_ = $dbalias{$_} || $_;
		if (m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$_ = uc($3) . ".sites." . $1 . '_' . $monthhash{$2};
			$monthhash{$2} or die "Error: the '$_' operation does not contain a valid month abbreviation\n";
		}
	}
	for (@db_names) {
		m/^mce(\d+way)$/ and $_ = "phastConsElements$1";		#for backward compatibility with older versions of ANNOVAR
	}
	return @db_names;
}

# Check whether the database file is present (and issue error message if not)
sub checkFileExistence {    
	my @db_names = @_;
	@genericdbfile1 = @genericdbfile;
	@vcfdbfile1 = @vcfdbfile;
	@bedfile1 = @bedfile;
	@gff3dbfile1 = @gff3dbfile;
	for my $i (0 .. @db_names-1) {
		my $dbfile;
		if ($db_names[$i] eq 'generic') {
			$dbfile = File::Spec->catfile ($dbloc, shift @genericdbfile1);
		} elsif ($db_names[$i] eq 'vcf') {
			$dbfile = File::Spec->catfile ($dbloc, shift @vcfdbfile1);
		} elsif ($db_names[$i] eq 'bed') {
			$dbfile = File::Spec->catfile ($dbloc, shift @bedfile1);
		} elsif ($db_names[$i] eq 'gff3') {
			$dbfile = File::Spec->catfile ($dbloc, shift @gff3dbfile1);
		} else {
			$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$db_names[$i].txt");
		}
		-f $dbfile or die "Error: the required database file $dbfile does not exist.\n";
	}
	@genericdbfile1 = @genericdbfile;
	@vcfdbfile1 = @vcfdbfile;
	@bedfile1 = @bedfile;
	@gff3dbfile1 = @gff3dbfile;
}

=head1 SYNOPSIS

 table_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --protocol <string>		comma-delimited string specifying database protocol
            --operation <string>	comma-delimited string specifying type of operation
            --outfile <string>		output file name prefix
            --buildver <string>		genome build version (default: hg18)
            --remove			remove all temporary files
            --(no)checkfile		check if database file exists (default: ON)
            --genericdbfile <files>	specify comma-delimited generic db files
            --gff3dbfile <files>	specify comma-delimited GFF3 files
            --bedfile <files>		specify comma-delimited BED files
            --vcfdbfile <files>		specify comma-delimited VCF files
            --otherinfo			print out otherinfo (infomration after fifth column in queryfile)
            --onetranscript		print out only one transcript for exonic variants (default: all transcripts)
            --nastring <string>		string to display when a score is not available (default: null)
            --csvout			generate comma-delimited CSV file (default: tab-delimited txt file)
            --argument <string>		comma-delimited strings as optional argument for each operation
            --tempdir <dir>		directory to store temporary files (default: --outfile)
            --vcfinput			specify that input is in VCF format and output will be in VCF format
            --dot2underline		change dot in field name to underline (eg, Func.refGene to Func_refGene)
            --thread <int>		specify the number of threads to be used in annotation
            --maxgenethread <int>	specify the maximum number of threads allowed in gene annotation (default: 6)
            

 Function: automatically run a pipeline on a list of variants and summarize
 their functional effects in a comma-delimited file, or to an annotated VCF file
 if the original input is a VCF file
 
 Example: table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,dbnsfp30a -operation g,r,f -nastring . -csvout
          table_annovar.pl example/ex2.vcf humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,dbnsfp30a -operation g,r,f,f,f,f,f,f,f -nastring . -vcfinput
                  
 Version: $Date: 2016-02-01 00:11:18 -0800 (Mon,  1 Feb 2016) $

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
(gene), r (region) or f (filter).

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--(no)checkfile>

the program will check if all required database files exist before execution of annotation

=item B<--genericdbfile>

specify the genericdb files used in -dbtype generic. Note that multiple comma-
delimited file names can be supplied.

=item B<--gff3dbfile>

specify the GFF3 dbfiles files used in -dbtype gff3. Note that multiple comma-
delimited file names can be supplied.

=item B<--bedfile>

specify the GFF3 dbfiles files used in -dbtype bed. Note that multiple comma-
delimited file names can be supplied.

=item B<--vcfdbfile>

specify the VCF dbfiles files used in -dbtype vcf. Note that multiple comma-
delimited file names can be supplied.

=item B<--otherinfo>

print out otherinfo in the output file. "otherinfo" refers to all the 
infomration after fifth column in the input queryfile.

=item B<--onetranscript>

print out only one random transcript for exonic variants. By default, all 
transcripts are printed in the output.

=item B<--nastring>

string to display when a score is not available. By default, empty string is 
printed in the output file.

=item B<--csvout>

generate comma-delimited CSV file. By default, tab-delimited text file is generated.

=item B<--argument>

a comma-separated list of arguments, to be supplied to each of the protocols. 
This list faciliates customized annotation procedure for each protocol.

=item B<--tempdir>

specify the directory location for storing temporary files used by 
table_annovar. This argument is especially useful in a cluster computing 
environment, so that temporary files are written to local disk of compute nodes, 
yet results files are written to possibly remote hosts.

=item B<--vcfinput>

specify that input is in VCF format and output will be in VCF format. if you 
want to generate a tab-delimited output or comma-delimited output file, you must 
use convert2annovar to generate an ANNOVAR input file first.

=item B<--dot2underline>

change dot in field name to underline (eg, Func.refGene to Func_refGene), which 
is useful for post-processing of the results in some software tools that cannot 
handle dot in field names.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

The table_annovar.pl program is designed to replace summarize_annovar.pl in 
earlier version of ANNOVAR. Basically, it takes an input file, and run a series 
of annotations on the input file, and generate a tab-delimited output file, 
where each column represent a specific type of annotation. Therefore, the new 
table_annovar.pl allows better customization for users who want to annotate 
specific columns.

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.

=cut
