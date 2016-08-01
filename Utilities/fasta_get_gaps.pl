#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie K
# version :  see below
# email   :  4urelie.k@gmail.com
# Purpose :  Read a genome in fasta formet to get gap coordinates, with output formatted as UCSC gap track
######################################################
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $scriptname = "fasta_get_gaps.pl";
my $version = "1.1";
my $changelog = "
#	- v1.0 = Apr 2013
#	- v1.1 = Jul 2016
#           Bug fix, lower cases n were not taken in account...
";


my $usage = "\nUSAGE (version $version):
    perl $scriptname -i <file.fasta> [-mlen <nt>] [-l] [-v] [-h]
	
    Will output a gap file (stretches on Ns), in the format described at: 
    http://genome.ucsc.edu/goldenPath/datorg.html
    columns are :
    <chromosome/ctg> <start-in-ctg> <end-in-ctg> <number> N <number-of-Ns> <kind> <bridged?>
	
    MANDATORY ARGUMENTS:
    -i,--in  (STRING) = fasta file
	
    OPTIONAL ARGUMENTS:
    -m,--mlen   (INT) = minimum length of a stretch of Ns to be considered as an assembly gap
	                    Any stretch of Ns > XX nt will be considered as a gap 
	                    [default = 100nt]
    -v,--v     (BOOL) = print version
    -c,--chlog (BOOL) = print change log (updates)
    -h,--help  (BOOL) = print this usage
\n";


#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($in,$v,$chlog,$help);
my $mlen = 100;
my $opt_success = GetOptions(
			 	  'in=s'      => \$in,
			 	  'm=s'       => \$mlen,
			 	  'chlog'     => \$chlog,
			 	  'version'   => \$v,
			 	  'help'      => \$help,);

#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if $v;
die $changelog if ($chlog);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$usage" if ($help || ! $in);
die "\n -f $in is not a fasta file?\n\n" unless ($in =~ /\.fa|fasta|fas$/);
die "\n -l $in does not exist?\n\n" if (! -e $in);
die "\n\tERROR: Minimum gap length (option -mlen) needs to be an integer\n$usage" if ($mlen !~ /\d/);


#-----------------------------------------------------------------------------
#------------------------------ NOW MAIN -------------------------------------
#-----------------------------------------------------------------------------

#index genome
my $fasta = Bio::SeqIO->new(-file => $in, -format => "fasta") or print LOG "ERROR: Failed to create SeqIO object from $in $!\n";

#get gaps
my $gaps = "$in.gaps.min$mlen.tab";
die "\n\tERROR: gap file already exists\n\n" if (-e $gaps);

my $log = "$in.gaps.min$mlen.log";
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
		$nt=uc($nt);
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
