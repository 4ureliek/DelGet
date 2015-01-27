#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie K
# date    :  v1.0, Apr 2013
# email   :  4urelie.k@gmail.com
# Purpose :  Read a genome in fasta formet to get gap coordinates, with output formatted as UCSC gap track
#			 Any stretch of Ns > XX nt will be considered as a gap [default = 100nt]
#
# 			 GAP FILE: http://genome.ucsc.edu/goldenPath/datorg.html
#					   it represents a gap it has the form:
#			           <chromosome/ctg> <start-in-ctg> <end-in-ctg> <number> N <number-of-Ns> <kind> <bridged?>
#
######################################################
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $ArgN = @ARGV;

my $usage = "\nUSAGE:
    perl scriptname.pl [-mlen <nt>] <genome.fasta>
	
	MANDATORY ARGUMENTS:
	<genome.fasta> = genome or any set of fasta sequences 
	
	OPTIONAL ARGUMENTS:
	-mlen <nt> : minimum length of a stretch of Ns to be considered as an assembly gap\n\n";

if ($ArgN < 1) {
	print $usage;
	exit;
}

my $Mlen;
my $mlen;
GetOptions ('mlen' => \$Mlen);
if ($Mlen){
	$mlen = shift @ARGV or die "$usage";
	die "\n\tERROR: Minimum gap length (option -mlen) needs to be an integer\n$usage" if ($mlen !~ /\d/);
} else {
	$mlen = 100;
}

#index genome
my $genome = shift @ARGV or die "$usage";
my $fasta = Bio::SeqIO->new(-file => $genome, -format => "fasta") or print LOG "ERROR: Failed to create SeqIO object from $genome $!\n";

#get gaps
my $gaps = "$genome.gaps.min$mlen.tab";
die "\n\tERROR: gap file already exists\n\n" if (-e $gaps);

my $log = "$genome.gaps.min$mlen.log";
open LOG, ">$log" or print LOG "ERROR: could not create file $log $!\n";
print LOG "Extracting gaps of:\n\n";

open GAPS, ">$gaps" or print LOG "ERROR: could not create file $gaps $!\n";
print GAPS "#bin\tchrom\tchromStart\tchromEnd\tix\tn\tsize\ttype\tbridge\n";
while( my $sequence = $fasta->next_seq() ) {
	my $head = $sequence->display_id;
	print LOG "$head\t";
	my $seq = $sequence->seq;
	my @seq = split("",$seq);
	my $Nstart = 0;
	my $Nlen = 0;
	NUCL: for (my $i = 0; $i <= $#seq; $i++){
		#go to next nt if this is not a gap to extend or first nt after a gap
		my $nt = $seq[$i];
		next NUCL if (($nt ne "N") && ($Nlen == 0));
		
		#loop on seq, if N => N start, then if following is N carry on until stops => Nend = Store in %gaps with ID = $Nstart
		if (($nt eq "N") && ($Nlen == 0)) { #1) first N in of the gap => get $Nstart + len = 1
			$Nstart = $i+1;
			$Nlen = 1;
		} elsif (($nt eq "N") && ($Nlen > 0)) { #2) at least second N => extend length
			$Nlen++;
		} elsif (($nt ne "N") && ($Nlen > 0)) { #3) len >0 => previous nt was a gap => print it and reinitialize values
			my $Nend = $Nlen + $Nstart - 1;
			print GAPS ".\t$head\t$Nstart\t$Nend\t.\t.\t$Nlen\t.\t.\n" unless ($Nlen < $mlen);
			$Nstart = 0;
			$Nlen = 0;
		}
	}
	print LOG "...done\n";
}
print LOG "\nGaps extracted, see $gaps\n";
close GAPS;
close LOG;

exit;

#output columns:
#bin	chrom	chromStart	chromEnd	ix	n	size	type		bridge
#585	GL429767	57033	57202		2	N	169		fragment	yes
