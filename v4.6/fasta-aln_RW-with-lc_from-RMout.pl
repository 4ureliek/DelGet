#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below / see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  DelGet util - Use an aligned file and a Repeat Masker outputs file (of the same sequences that were aligned) to rewrite sequences masked in lowercases
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::SeqIO;

my $version = "2.0";
my $scriptname = "fasta-aln_RW-with-lc_from-RMout.pl";
my $changelog = "
#	- v1.0 = May 2013
#	- v2.0 = 08-09 Sept 2016
#            Complete re writing - subroutines + simplification
#            Use the .masked instead of the .out
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname -i <directory> [-v] [-h] [-l]
    			
    PURPOSE:
    Primarily a DelGet util (https://github.com/4ureliek/DelGet), but could be more general;
    The script detects repeat masker output files .fa.masked
    and will rewrite the masked regions in lowercases for the .fa.align.fa files
    Note that it expects masked sequences as X
		
    MANDATORY ARGUMENTS:	
     -i,--in      => (STRING) directory with files to be processed
                              (requries the repeat masker output .masked and the aligned .fa.align.fa file)
 
    OPTIONAL ARGUMENTS
     -v,--v       => (BOOL)   verbose, makes the script talk to you
     -v,--v       => (BOOL)   print version if only option
     -h,--help    => (BOOL)   print this usage
     -l,--log     => (BOOL)   print the change log (updates)
\n";


################################################################################
### Get arguments/options
### check some of them, print details of the run if verbose chosen
################################################################################
my ($in,$help,$chlog,$v);
GetOptions ('in=s'      => \$in,  
            'log'       => \$chlog, 
            'help'      => \$help, 
            'v'         => \$v);

#check step for options
die "\n $scriptname version: $version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n Please provide input directory (-i, use -h to see the usage)\n\n" if (! $in);
die "\n directory -i $in does not exist?\n\n" if (($in !~ /,/) && (! -e $in));
$in = $1 if ($in =~ /^(.*)\/$/); #remove the / at the end if any

################################################################################
### MAIN
################################################################################
print STDERR "--------------------------------------------------\n" if ($v);
print STDERR " --- Obtaining list of files to deal with\n" if ($v);
my $masked = input_files($in,$v);

########################################
# Loop on these regions
########################################
foreach my $reg (@{$masked}) {	
	my $aln = "$reg.fa.align.fa";
	print STDERR " --- Dealing with $reg ($aln)\n" if ($v);
	
	my $masked = "$reg.fa.masked";
	print STDERR "     Getting masked positions from $masked\n" if ($v);
	my $rm = get_masked($masked,"X",$v);
	
	my $alnrm = "$reg.fa.align.masked.fa";
	print STDERR "     Rewriting sequences with masked positions as lowercases in $alnrm\n" if ($v);	
	aln_rw_lc($aln,$rm,$alnrm,$v);
}

print STDERR " --- Done, see output files *.masked.fa in $in\n" if ($v);
exit;





################################################################################
### SUBROUTINES
################################################################################
#----------------------------------------------------------------------------
# get list of input files
# my $masked = input_files($in,$v);
#----------------------------------------------------------------------------
sub input_files {
	my($in,$v) = @_;
	my @masked = `ls $in/*masked`;
	my @reg = ();
	foreach my $reg (@masked) {
		chomp($reg);
		$reg =~ s/\.fa\.masked//;
		push(@reg,$reg);
	}
	return (\@reg);
}

#----------------------------------------------------------------------------
# get masked positions
# my $rm = get_masked($masked,"X",$v);
#----------------------------------------------------------------------------
sub get_masked {
	my($masked,$t,$v) = @_;
	my %posi = ();
	my $aln_in = Bio::SeqIO->new(-format => 'Fasta', -file=>$masked) or confess "     \nERROR (sub get_masked): could not create Bio::SeqIO object from $masked $!\n";
	while( my $seq = $aln_in->next_seq() ) {
		my $id = $seq->display_id;
		my $alnlen = $seq->length();
		my @seq = split(//,$seq->seq);
		for (my $i=0; $i<$alnlen; $i++) {
			($seq[$i] eq $t)?($posi{$id}[$i]="y"):($posi{$id}[$i]="n");
		}
	}
	unlink $aln_in;
	return (\%posi);
}

#----------------------------------------------------------------------------
# rewrite alignment with masked nt in lower case
# aln_rw_lc($aln,$rm,$alnrm,$v);
#----------------------------------------------------------------------------
sub aln_rw_lc {
	my($aln,$rm,$alnrm,$v) = @_;
	my $aln_io = Bio::SeqIO->new(-format => 'Fasta', -file=>$aln) or confess "     \nERROR (sub aln_rw_lc): could not create Bio::SeqIO object from $aln $!\n";
	open(my $fho, ">", $alnrm) or confess "     \nERROR (sub aln_rw_lc): could not open to write $alnrm $!\n";	
	my $t = "X"; #because I masked with X - this script won't work otherwise
	SEQ: while( my $seq = $aln_io->next_seq() ) {			
		my $id = $seq->display_id;
		print STDERR "      -> $id\n" if ($v);
		print $fho ">$id\n";
		my $alnseq = $seq->seq;	
		#print the sequence if nothing was masked and go to the next
		print $fho "$alnseq\n" if (! $rm->{$id});
		next SEQ if (! $rm->{$id});
		#now loop on the sequence to rewrite accordingly
		my $len = $seq->length();	
		my $c = 0;
		my @seq = split(//,$alnseq);
		my $masked = ""; #to store the sequence and limit the IO writings
		for (my $i=0;$i<$len;$i++) {
			if ($seq[$i] eq "-") {
				$masked=$masked.$seq[$i];
			} else {
				#not a gap; $c did not get incremented while gaps => real coordinate in $rm
				($rm->{$id}[$c] eq "y")?($masked=$masked.lc($seq[$i])):($masked=$masked.$seq[$i]);
				$c++;
			}
		}
		print $fho "$masked\n";
	}
	close $fho;
	return 1;
}




# 
# 
# #----------------------------------------------------------------------------
# # DEPRECATED: get intervals
# # my $gaps = fasta_get_intervals($aln,"-",$v); <= deprecated, no need
# # my $rm = fasta_get_intervals($masked,"X",$v);
# #----------------------------------------------------------------------------
# sub fasta_get_intervals {
# 	my($aln,$t,$v) = @_;
# 	my %posi = ();
# 	my $aln_in = Bio::SeqIO->new(-format => 'Fasta', -file=>$aln) or confess "     \nERROR (sub get_gaps): could not create Bio::SeqIO object from $aln $!\n";
# 	while( my $seq = $aln_in->next_seq() ) {
# 		my $header = $seq->display_id;
# 		my $alnlen = $seq->length();
# 		my @seq = split(//,$seq->seq);
# 		my $start;
# 		for (my $n=1;$n<$alnlen;$n++) {
# 			#gap opening: get start
# 			$start = $n if (($n == 1) && ($seq[$n-1] eq $t)); #just case of first nt of aln
# 			$start = $n+1 if (($seq[$n] eq $t) && ($seq[$n-1] ne $t));
# 			if (($seq[$n-1] eq $t) && ($seq[$n] ne $t)) { 
# 			#gap ending: get end + print coords
# 				my $end = $n+1;
# 				my $len = $end - $start;
# 				$posi{$header}.="$start\t$len\t$end\t";
# 			}
# 		}
# 	}
# 	unlink $aln_in;
# 	return (\%posi);
# }
# 
# 



