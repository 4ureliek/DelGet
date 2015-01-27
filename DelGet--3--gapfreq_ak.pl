#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie K
# version :  1.4 (see below)
# email   :  4urelie.k@gmail.com
# Pupose  :  For all files *.align.fa in provided path, get gap counts + length
#				=> compare gap profiles
#				Note that it has be written for 3 species aln; theoretically works for more (loops will loop) but not tested
#				However REQUIREMENT is that all aln have the same number of sequences (assumed 0..n being always same order)
######################################################
# UPDATES
#	- v1.0 = 07 May 2013
#	- v1.1 = 09 Oct 2013
#		what's printed + output
#	- v1.2 = 15 Oct 2013
#		Corrected missing full path to open files (integration into pipeline)
#		Corrected some var names / values
#	- v1.3 = 18 Oct 2013
#		"align" instead of "muscle" in names
#	- v1.4 = 07 Nov 2013
#		Do not consider that the A at : TGC------A------TGC is a gap ending, but go to next. 
#		Coordinates of start and end not changed, but length of the gap is corrected
##########################################################
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Array::Transpose::Ragged qw/transpose_ragged/;
my $version = "v1.4";

my $usage = "Usage:
	perl <scriptname.pl> <path_to_all> <spID1,speID2,spID3>

	With:
	<path_to_all> containing all files to be checked, with names as *.align.fa
	<spID1,speID2,spID3> being the order of species IDs (that are in seq names in alignment, as spID1_gi|...|etc), as:
		      |-----sp2
		|-----|
		|     |-----sp3
		|  
		|-----------sp1 [first ID needs to be the outgroup one]
		
	Typically:
	perl script.pl /data/user/del_analysis/aln_seq Efus,Mluc,Mdav\n\n";
	
my $path = shift @ARGV or die "$usage";
my $spIDs = shift @ARGV or die "$usage";
my @spIDs = split(/,/,$spIDs);

print "Script gapfreq version $version started:\n";
print "@spIDs\n";

##########################################################################################################
# Get list of files
##########################################################################################################
my @List = `ls $path`;
my @Files;
foreach my $file (@List) {
	chomp $file;
	my $full_file = "$path/$file";
	push(@Files,$full_file) if ($file =~ m/^.*\.align\.fa$/);
}

##########################################################################################################
# Get gaps for each aln => one file per sequence/specie in each alignement
##########################################################################################################
my %totlen = ();
foreach my $file (@Files) {
	print "--- $file in progress\n";
	#####################################################
	# Extract sequences of align aln files => in matrix + prep GAPOUT
	#####################################################
	print "    Getting sequences from alignement\n";
	my %seqs = ();
	my $length;
	my $aln = Bio::SeqIO->new(-file => $file, -format => "fasta") or die "Failed to create SeqIO object from $file $!\n";
	while( my $seq = $aln->next_seq() ) {
		my $ID = $seq->display_id;
		$ID =~ s/^([^_][a-zA-Z1-9]+)_.*$/$1/;
		my @sequence = split(//,$seq->seq);
		$seqs{$ID} = \@sequence;
		$length = $seq->length;
		$totlen{$ID}+=$length;
	}
	
	
	#####################################################
	# Get gaps and print their coordinates (all files in same gap file => treat all gaps at the same time)
	#####################################################
	print "    Getting gap coordinates in files\n";
	my %start = ();
	my %gap_nb = ();
	my %skipped = ();
	for (my $n=1;$n<($length-1);$n++) {
		for (my $i=0;$i<=$#spIDs;$i++) {
			#nt for each species at each position is ($seqs{$spIDs[$i]}->[$n]
			#print in a file gap start-end couples, one per line, one file per species [to keep track]
			my $gaps_out = "$path/_$spIDs[$i].gaps.bed";
			open(GAPOUT,">>$gaps_out") or die "\t    ERROR - can not open file $gaps_out $!";
			my $region = $file;
			$region =~ s/^reg(.*)\.fa\.align\.fa$/reg$1/;
			if (($seqs{$spIDs[$i]}->[$n] eq "-") && ($seqs{$spIDs[$i]}->[$n-1] ne "-") && ($seqs{$spIDs[$i]}->[$n-2] ne "-")) { 
			# "-" at position $n => this is a gap opening; IF $n-1 ne "-" AND IF $n-2 ne "-" [ie situation TGC-----A-----TGC]
				$start{$spIDs[$i]} = $n+1;
			}
			if (($seqs{$spIDs[$i]}->[$n-1] eq "-") && ($seqs{$spIDs[$i]}->[$n] ne "-") && ($start{$spIDs[$i]}) && ($seqs{$spIDs[$i]}->[$n+1] eq "-")) {
			#this is to correct length in case there is a case like A in TGC-----A-----TGC, because won't be considered as a gap ending.
				($skipped{$spIDs[$i]})?($skipped{$spIDs[$i]}+=1):($skipped{$spIDs[$i]}=1); #remember how many
			}			
			if (($seqs{$spIDs[$i]}->[$n-1] eq "-") && ($seqs{$spIDs[$i]}->[$n] ne "-") && ($start{$spIDs[$i]}) && ($seqs{$spIDs[$i]}->[$n+1] ne "-")) { 
			#this is a gap ending, unless gap is on begining of alignement - I don't want to count these since I don't know their length
			#and unless it's a case like A in TGC-----A-----TGC (if one or more was seen, start coordinate was "corrected")
				my $end = $n;
				my $len;
				if ($skipped{$spIDs[$i]}) {
					$len = $end - $start{$spIDs[$i]} +1 - $skipped{$spIDs[$i]};
					$skipped{$spIDs[$i]} = 0; #need to be reinitialized to 0 for next gap
				} else {
					$len = $end - $start{$spIDs[$i]} +1;
				} 
				$gap_nb{$spIDs[$i]}++;
				print GAPOUT "$region\t$start{$spIDs[$i]}\t$end\t$gap_nb{$spIDs[$i]}\t.\t+\t$len\t$length\n";
				
			}
			close GAPOUT;
		}
	}
}
print "--- done\n";


##########################################################################################################
# Now analyse gap files
##########################################################################################################
# Get deletions specific to species 1 or 2
##########################################################################################################
# 1) => bedtools on gap files
# => get only gaps that are not in the 2 others
# => but allow a certain flexibility though (ie should be considered independant gaps => specie spe if they are not same boundaries
#    even though of course not possible to know if there was 1 first event shared and then another event in only one specie
#	 consider same gap if overlap is >85%
#	 this min overlap is for BOTH sense, ie if BIG gap in one specie, and a small one in the other => indep
#	 still needs to subtract the other species' stuff, if TE insterted (or there) at deletion point. Note that TE could have been shared in the first place and deleted for real, but I prefer underestimate
print "--- parsing gaps\n";

my $outgroup = "$path/_$spIDs[0].gaps.bed";
my $gaps_sp1 = "$path/_$spIDs[1].gaps.bed";
my $gaps_sp2 = "$path/_$spIDs[2].gaps.bed";

my $sub_1out = "$path/_$spIDs[1].gaps.specific.bed";
system "intersectBed -a $gaps_sp1 -b $gaps_sp2 -f 0.85 -r -v | intersectBed -a stdin -b $outgroup -f 0.85 -r -v | subtractBed -a stdin -b $gaps_sp2 -f 0.1 | subtractBed -a stdin -b $outgroup -f 0.1 > $sub_1out";

my $sub_2out = "$path/_$spIDs[2].gaps.specific.bed";
system "intersectBed -a $gaps_sp2 -b $gaps_sp1 -f 0.85 -r -v | intersectBed -a stdin -b $outgroup -f 0.85 -r -v | subtractBed -a stdin -b $gaps_sp1 -f 0.1 | subtractBed -a stdin -b $outgroup -f 0.1 > $sub_2out";
my @spegapfiles = ($sub_1out,$sub_2out);



# 2) => Now parse subtracted file => separate small and big + calculate amounts, print distribution
foreach my $gapfile (@spegapfiles) {
	my $ID = $gapfile;
	$ID =~ s/^$path\/_([a-zA-z1-9]+)\.gaps\.specific\.bed$/$1/;
	#create ouput files
	my $spesmall = "$gapfile.small.tab";
	open(SPESMALL,">$spesmall") or die "\t    ERROR - can not open file $spesmall $!"; 
	my $spebig = "$gapfile.big.tab";
	open(SPEBIG,">$spebig") or die "\t    ERROR - can not open file $spebig $!"; 
	#loop on gap file
	open(SPE,"<$gapfile") or die "\t    ERROR - can not open file $gapfile $!";
	my %smallgapcount = ();
	my $smallgapslen = 0;
	my %biggergapcount = ();
	my $biggapslen = 0;
	while (<SPE>) {
		#0		1	2	3	4	5	6	7
		#reg2-4	827	828	2	.	+	1	11296
		chomp (my $line = $_);
		my @line = split(/\t/,$line);
		my $IDfull = $line[0]."_".$line[3];
		my $gaplen = $line[2] - $line[1] + 1;
	
		#separate + count events (gap ID)
		if ($gaplen <=  30) {
			$smallgapcount{$IDfull}++;
			$smallgapslen+= $gaplen;
			print SPESMALL "$line\n";
		} else {
			$biggergapcount{$IDfull}++;
			$biggapslen+= $gaplen;
			print SPEBIG "$line\n";
		}
	}
	close SPE;
	close SPESMALL;
	close SPEBIG;
	my $smallgaps = 0;
	foreach my $small (sort keys %smallgapcount)  {
		$smallgaps++;
	}
	my $biggergaps = 0;
	foreach my $big (sort keys %biggergapcount)  {
		$biggergaps++;
	}

	print "    $gapfile:\n";
	print "    Total length of alignements\t$totlen{$ID}\n";
	print "    - Total number of spe gaps 1-30nt\t$smallgaps\n";
	print "    - Total length of spe gaps 1-30nt\t$smallgapslen\n";
	print "    - Total number of spe gaps >30nt\t$biggergaps\n";
	print "    - Total length of spe gaps >30nt\t$biggapslen\n";
}
print "--- done\n";
##########################################################################################################
exit;

















