#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $REVISION = '$Revision: 331c86ef86c8906f67007b0cc2a2c64f93532b4c $';
our $DATE =	'$Date: 2015-12-14 18:18:18 -0800 (Mon, 14 Dec 2015) $';  
our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org> $';

our ($verbose, $help, $man, $includesnp, $mrnaseq, $onlyAltering);

our ($evffile, $genefile, $fastafile);

our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'includesnp'=>\$includesnp, 'mrnaseq'=>\$mrnaseq, 'onlyAltering'=>\$onlyAltering) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 3 or pod2usage ("Syntax error");

($evffile, $genefile, $fastafile) = @ARGV;

open (EVF, $evffile) or die "Error: cannot read from evffile $evffile: $!\n";
open (GENE, $genefile) or die "Error: cannot read from genefile $genefile: $!\n";
open (FASTA, $fastafile) or die "Error: cannot read from fastafile $fastafile: $!\n";

my (@queue, %need_trans);
while (<EVF>) {
	s/[\r\n]+$//;
	m/^line\d+/ or die "Error: invalid record found in exonic_variant_function file $evffile (line number expected): <$_>\n";
	my @field = split (/\t/, $_);
	$field[2] =~ m/^\w+:(\w+):wholegene/ and next;
	$field[2] =~ m/unknown/i and next;
	
	#ensembl gene name may contain - and ., where as UCSC transcript name may contain '.'.
	$field[2] =~ m/^[\w\-\.]+?:([\w\.]+?):exon\d+:c.([\w>]+)/ or warn "Error: invalid record found in exonic_variant_function file (exonic format error): <$_>\n" and next;
	my ($transcript, $cchange) = ($1, $2);
	my ($start, $end, $ref, $obs);
	if ($cchange =~ m/^([ACGTacgt])(\d+)([ACGTacgt])$/) {
		if ($includesnp) {
			($start, $end, $ref, $obs) = ($2, $2, $1, $3);
		} else {
			next;
		}
	} elsif ($cchange =~ m/^(\d+)([ACGTacgt])>([ACGTacgt])$/) {
		if ($includesnp) {
			($start, $end, $ref, $obs) = ($1, $1, $2, $3);
		} else {
			next;
		}
	} elsif ($cchange =~ m/^(\d+)_(\d+)delins(\w+)/) {	#block substitution
		($start, $end, $obs) = ($1, $2, $3);
	} elsif ($cchange =~ m/^(\d+)del\w+/) {		#single base deletion
		($start, $end, $obs) = ($1, $1, '');
	} elsif ($cchange =~ m/^(\d+)_(\d+)del(\w*)/) {	#multi-base deletion
		($start, $end, $obs) = ($1, $2, '');
	} elsif ($cchange =~ m/^(\d+)_(\d+)ins(\w+)/) {	#insertion
		($start, $end, $ref, $obs) = ($1, $1, 0, $3);		#if end is equal to start, this is an insertion
	} elsif ($cchange =~ m/^(\d+)dup(\w+)/) {	#insertion
		($start, $end, $ref, $obs) = ($1, $1, 0, $2);
	} elsif ($cchange =~ m/^(\d+)_(\d+)(\w+)/) {	#non-frameshift substitution
		($start, $end, $obs) = ($1, $2, $3);
	} else {
		die "Error: invalid coding change format: <$cchange> within <$_>\n";
	}
	push @queue, [$field[0], $transcript, $start, $end, $ref, $obs, $cchange];
	$need_trans{$transcript}++;
}


my (%mrnastart, %mrnaend);
while (<GENE>) {
	s/[\r\n]+$//;
	my @field = split (/\t/, $_);
	@field >= 11 or die "Error: invalid record found in gene file (>=11 fields expected): <$_>\n";
	$field[0] =~ m/^\d+$/ and shift @field;		#refGene and ensGene has bin as the first column

	my ($name, $strand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend) = @field[0, 2, 3, 4, 5, 6, 8, 9];
	$need_trans{$name} or next;
	
	my ($mrnastart, $mrnaend);
	
	#next we need to make sure that there is no intron between transcription start and translation start (this is rare but it happens when cdsstart is not in the first exon)
	my @exonstart = split (/,/, $exonstart);
	my @exonend = split (/,/, $exonend);
	
	$txstart++;
	$cdsstart++;
	@exonstart = map {$_+1} @exonstart;

	if ($strand eq '+') {
		#<---->-----<--->----<------>----<----->-----<--->
		#             **********
		my $intron = 0;
		for my $i (0 .. @exonstart-1) {
			$i and $intron += ($exonstart[$i]-$exonend[$i-1]-1);
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				$mrnastart = $cdsstart-$txstart+1-$intron;
			}
			if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
				$mrnaend = $cdsend-$txstart+1-$intron;
			}
			
		}
	} elsif ($strand eq '-') {
		#<---->-----<--->----<------>----<----->-----<--->
		#             **********
		my $intron = 0;
		for (my $i=@exonstart-1; $i>=0; $i--) {
			$i<@exonstart-1 and $intron += ($exonstart[$i+1]-$exonend[$i]-1);
			if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
				$mrnastart = $txend-$cdsend+1-$intron;
			}
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				$mrnaend = $txend-$cdsstart+1-$intron;
			}
			
		}
	}
			

	$mrnastart{$name} = $mrnastart;
	$mrnaend{$name} = $mrnaend;
}

my (%mrnaseq);
my ($curname, $curseq);

while (<FASTA>) {
	s/[\r\n]+$//;
	if (m/^>([\w\.]+)/) {
		if ($curseq) {
			$mrnaseq{$curname} = $curseq;
		}
		$curname = $1;
		$curseq = '';
	} else {
		$curseq .= $_;
	}
	$curseq and $mrnaseq{$curname} = $curseq;	#process the last sequence
}

#process each element in the queue
for my $i (0 .. @queue-1) {
	my ($line, $transcript, $start, $end, $ref, $obs, $cchange) = @{$queue[$i]};
	$verbose and print STDERR "NOTICE: Processing $line with $transcript, $start, $end, $obs, $cchange\n";
	if (not defined $mrnaseq{$transcript}) {
		print STDERR "WARNING: cannot find mRNA sequence for $transcript in the fastafile $fastafile\n";
		next;
	}
	if (not defined $mrnastart{$transcript}) {
		print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription start site\n";
		next;
	}
	if (not defined $mrnaend{$transcript}) {
		print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription end site\n";
		next;
	}
	
	if ($end > length ($mrnaseq{$transcript})) {
		print STDERR "ERROR: transcript end ($mrnaend{$transcript}) for $transcript is longer than transcript length ${\(length ($mrnaseq{$transcript}))}, skipping this transcript\n";
		next;
	}
	
	
	if ($mrnaseq) {
		my $utr5 = substr ($mrnaseq{$transcript}, 0, $mrnastart{$transcript}-1);
		my $utr3 = substr ($mrnaseq{$transcript}, $mrnaend{$transcript});
		my ($mrna1, $mrna2);
		$mrna1 = $mrnaseq{$transcript};
		
		my $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}-1, $mrnaend{$transcript}-$mrnastart{$transcript}+1);
		my @dna = split (//, $dna);
		my ($protein1, $protein2);
		my $warning = '';
	
		if ($end > @dna) {
			print STDERR "ERROR in $line: end position of variant ($end) in $transcript is longer than coding portion length ${\(scalar @dna)}, skipping this transcript\n";
			next;
		}
	
		$protein1 = translateDNA ($dna);
		$protein1 =~ m/\*\w/ and $warning = '(WARNING: Potential FASTA sequence error!!!)';
		if ($start == $end and not $ref and $obs) {		#this is an insertion
			splice (@dna, $start, 0, $obs);
		} else {				#this is a substitution
			splice (@dna, $start-1, $end-$start+1, $obs);
		}
		$dna = join ('', @dna);
		
		$mrna2 = $utr5 . $dna . $utr3;
				
		$mrna1 =~ s/(.{100})/$1\n/g;
		print ">$line $transcript WILDTYPE $warning\n";
		print "$mrna1\n";
		$mrna2 =~ s/(.{100})/$1\n/g;
		print ">$line $transcript c.$cchange\n";
		print "$mrna2\n";
		
	} else {
	
		my $dna = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}-1, $mrnaend{$transcript}-$mrnastart{$transcript}+1);
		my @dna = split (//, $dna);
		my ($protein1, $protein2);
		my $warning = '';
	
		if ($end > @dna) {
			print STDERR "ERROR in $line: end position of variant ($end) in $transcript is longer than coding portion length ${\(scalar @dna)}, skipping this transcript\n";
			next;
		}
	
		$verbose and print STDERR "NOTICE: wild-type DNA sequence is $dna\n";
		$protein1 = translateDNA ($dna);
		$protein1 =~ m/\*\w/ and $warning = '(WARNING: Potential FASTA sequence error!!!)';
		if ($start == $end and not $ref and $obs) {		#this is an insertion
			splice (@dna, $start, 0, $obs);
			$start++;
			$end++;		#after insertion of $obs, the position of start/end increase by one
		} else {				#this is a substitution
			splice (@dna, $start-1, $end-$start+1, $obs);
		}
		
		my $aastart = int(($start-1)/3)+1;
		my $aaend1 = int (($end-1)/3)+1;
		my $aaend2 = int (($start+length($obs)-1-1)/3)+1;
		my ($function, $aachange);
		$aachange = '';		#assign default value, since synonymous variants do not have aachange
		
		$dna = join ('', @dna);
		$protein2 = translateDNA ($dna);
		$verbose and print STDERR "NOTICE:   mutated DNA sequence is $dna\n";
		
		if (substr ($protein2, $aastart-1, $aaend2-$aastart+1) =~ m/\*/) {
			$function = 'immediate-stopgain';
		} elsif ($end>=$mrnaend{$transcript}-2) {
			$function = 'immediate-stoploss';
		} else {
			if ($protein1 eq $protein2) {
				$function = 'synonymous';
				$onlyAltering and next;
			} else {
				$function = 'protein-altering';		#change nonsynonymous to protein-altering to reduce confusion among users 20150604
			}
		}
		
		#I changed the paragraph below 20150604 to reduce confusion among users: for insertions, I now identify the exact position that has change (rather thant the position where nucleotide changes)
		#$aachange = " (position $aastart-$aaend1 changed from " . substr ($protein1, $aastart-1, $aaend1-$aastart+1) . ' to ' . substr ($protein2, $aastart-1, $aaend2-$aastart+1). ')';
		for my $i (0 .. length($protein1)-1) {
			if (substr($protein1, $i, 1) ne substr($protein2, $i, 1)) {
				$aachange = " (position ${\($i+1)} changed from " . substr($protein1, $i, 1) . ' to ' . substr($protein2, $i, 1) . ")";
				last;
			}
		}
		
		
		$protein2 =~ s/\*.+/\*/;		#delete anything after the STOP codon
	
	
		$protein1 =~ s/(.{100})/$1\n/g;
		print ">$line $transcript WILDTYPE $warning\n";
		print "$protein1\n";
		$protein2 =~ s/(.{100})/$1\n/g;
		print ">$line $transcript c.$cchange $function $aachange\n";
		print "$protein2\n";
	}
}





sub translateDNA {
	my ($seq) = @_;
	my ($nt3, $protein);
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codon1{$1} or print "WARNING: invalid triplets found in DNA sequence to be translated: <$1> in <$seq>\n" and die;
		$protein .= $codon1{$1};
	}
	return $protein;
}



=head1 SYNOPSIS

 coding_change.pl [arguments] <exonic-variant-function-file> <gene-def-file> <fasta-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --includesnp		process SNPs (default is to process indels only)
            --mrnaseq			print out mRNA sequence (including UTR and coding)
            --onlyAltering		ignore synonymous SNPs (when -includesnp is on)

 Function: infer the translated protein sequence (or mRNA sequence) for exonic 
 frameshift mutations (or SNPs) identified by ANNOVAR
 
 Example: coding_change.pl ex1.avinput.exonic_variant_function humandb/hg19_refGene.txt humandb/hg19_refGeneMrna.fa
          coding_change.pl ex1.avinput.exonic_variant_function humandb/hg19_refGene.txt humandb/hg19_refGeneMrna.fa -includesnp -onlyAltering
 
 Version: $Date: 2015-12-14 18:18:18 -0800 (Mon, 14 Dec 2015) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--includesnp>

include sequences affected by SNPs in the output. By default, only sequences 
affected by exonic indels are printed out.

=item B<--mrnaseq>

print out mRNA sequences rather than protein sequences.

=item B<--onlyAltering>

when -includesnp is set, do not process synonymous variants and only print out
protein-coding alterations

=back

=head1 DESCRIPTION

This program will infer the protein sequence for frameshift mutations identified 
by ANNOVAR. Typically, for non-synonymous mutations, the annotate_variation.pl 
program in ANNOVAR will report the amino acid change at the position; however, 
for frameshift mutations which may affect a long stretch of amino acid, ANNOVAR 
will only give a frameshift annotation without printing out the new protein 
sequence. This program can take ANNOVAR exonic_variant_function file and attempt 
to infer the new protein sequence.

The <exonic-variant-function-file> is the one generated by
annotate_variation.pl. The <gene-def-file> is the gene definition file that is
used by annotate_variation.pl (which can be examined by looking at the LOG
message. The <fasta-file> is also used by annotate_variation.pl (which can be
examined by looking at the LOG message

=over 8

=item * B<Known Bug>

This program does not handle mitochondria mutations correctly yet.

=back

For questions, comments or bug reports, please contact me at 
$Author: Kai Wang <kai@openbioinformatics.org> $.

=cut