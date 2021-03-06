#!/usr/bin/perl -w

##########################################################
# Author  :  Aurelie Kapusta
# version :  3.0 (see below)
# email   :  4urelie.k@gmail.com
# PURPOSE :  General question at the base of writing this script = assess medium size deletions between species.
#			 See $usage and CONFIG FILE $help for details
##########################################################
# UPDATES
#	- v1.0 = 02 May 2013
#	- v2.0 = 14 May 2013
#		lots of bug fixing beetween v1.0 and v2.0
#		now part of a pipeline
#	- v2.1 = 20 May 2013
#		check that distance between anchors >=0
#	- v2.2 = 26 Jun 2013
#		possibility of loading previous output file, to avoid overlaps => append new regions to a previous run
#	- v2.3 = 19 Jul 2013
#		maximum anchor length = 150% input (and not fixed at 120nt)
#	- v2.4 = 08 Sep 2013
#		Corrected error at the step "#3) avoid overlap with any already picked region" when positions are randomized that caused regions to overlap
#				if (($three_end < ($posi + 2*$anch_len + $anch_dist - 3)) && ($three_end > ($posi - 2*$anch_len + $anch_dist + 3)))
#			is now:
#				if (($three_end < ($posi + 2*$anch_len + $anch_dist - 3)) && ($three_end > ($posi - 2*$anch_len - $anch_dist + 3)))
#	- v2.5 = 15 Oct 2013
#		change way of getting previous output
#		re integration into pipeline, regions obtained $randnb by $randnb => big loop with loading previous outputs etc moved into pipeline script
#		changes in $OKregions loading
#	- v2.6 = 06 Nov 2013
#		Regions were still overlapping despite loading $OKregions => %alreadyrand_full was not declared properly
#		Added version follow up
#	- v2.7 = 14-20 Nov 2013
#		Regions were still overlapping -> debugged and fixed it [for real this time]; 
#			(issue was not previous runs loading but during the run itself. Wrong value was stored (not end of region))
#		+ during debugging I ended up changing ways of checking for overlap and double checked that assembly gap overlap was not happening
#	- v2.8 = 07 Jan 2014
#		Try to avoid running without printing. 
#			=> All exists tests -> define tests
#			=> store of previous positions to avoid in a doc that will be close and reopened at each round
#	- v3.0 = 05 Feb 2014
#		Changes to use subroutines => lots of changes in the code
#	- v3.1 = 04 Mar 2014
#		Bug corrected for when Deletions.2 => was trying to open Deletions.0/_OKregions.all.tab instead of in Deletions.1
#   - v3.2 = 27 Jan 2015
#       Changes for Github first upload (@INC stuff)
#   - v3.3 = 13 May 2015
#       Few changes of no consequences while debugging
#       Set counter $failed_anchor_ref->{$r} to 0
#		in DelGet.pm, sub get_all_len_filtered 
#          => if ratio = 0, put it at 1 - otherwise it can never go through when many sequences and few are asked
#   - v3.4 = 26-27 July 2016
#       Bug fix in loading gaps (suroutine get_gaps, only last genome was loaded
#       Also the -- for the ratio was before filters => was going down before a succesfull region was picked
##########################################################
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN/Lib");
}
use Array::Unique;
use DelGet;
my $version = "v3.4";

####################################################################################################################
# Usage and load config file
####################################################################################################################
my $usage = "\nUSAGE:
	perl script.pl <config_file> <path>

	Typically:
		perl DelGet--1--get_regions.pl DelGet--1--CONFIG_manipname.pl Deletions.X
		
	This script is part of a pipeline, type \"perl DelGet.pl -help\" for more details\n\n";
	
my $config_file = shift @ARGV or die "$usage $!\n" ;
my $pathtemp = shift @ARGV or die "$usage $!\n" ;

#################################################################
# Variables and names
#################################################################
# Load config file
my $cf = DelGet::load_config($config_file);
#temporary fix until I replace the variable values in this script
my %cf = %{$cf};
my ($a_max_len, $randtot, $randnb, $a_per_round, $anch_dist, $anch_len, $maxratio, $minlen, $mingaplen, $multip, $BLATSOFT, $OKregions, $if_OKregions, $ALNSOFT, $RMSOFT, $path, $genone, $gentwo, $genthree, $genone_a, $gentwo_a, $genthree_a, $IDgen1, $IDgen2, $IDgen3, $gapone, $gaptwo, $gapthree, $mask, $masking_folders) = ($cf{'a_max_len'}, $cf{'randtot'}, $cf{'randnb'}, $cf{'a_per_round'}, $cf{'anch_dist'}, $cf{'anch_len'}, $cf{'maxratio'}, $cf{'minlen'}, $cf{'mingaplen'}, $cf{'multip'}, $cf{'BLATSOFT'}, $cf{'OKregions'}, $cf{'if_OKregions'}, $cf{'ALNSOFT'}, $cf{'RMSOFT'}, $cf{'path'}, $cf{'genone'}, $cf{'gentwo'}, $cf{'genthree'}, $cf{'genone_a'}, $cf{'gentwo_a'}, $cf{'genthree_a'}, $cf{'IDgen1'}, $cf{'IDgen2'}, $cf{'IDgen3'}, $cf{'gapone'}, $cf{'gaptwo'}, $cf{'gapthree'}, $cf{'mask'}, $cf{'masking_folders'});

# Get lists of files
my @genfiles = ($gentwo,$genthree);
my @genIDs = ("$IDgen2","$IDgen3");
my @gapfiles = ($gapone,$gaptwo,$gapthree);

#Folder for outputs
my $c=$pathtemp;
$c =~ s/^.*Deletions\.([0-9]+)$/$1/;

#$logone_fh FILE
my $log = "$pathtemp/_DelGet--1--get_regions.log";
open(my $logone_fh, ">",$log) or die "\t    ERROR - can not open to write log file $log $!\n";
print $logone_fh "\nScript DelGet--1--get_regions.pl version $version started:\n";
print $logone_fh "\n---GET REGIONS----------------------------------------------\n    Run number $c\n";

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

###########################################################
# DEAL WITH GAP FILES
###########################################################
print $logone_fh "\n---GAP FILES----------------------------------------------\n";
#0		1							2			3			4	5	6
#bin	chrom						chromStart	chromEnd	ix	n	size	type	bridge
#.		gi|432120764|gb|KB097841.1|	268			451			.	.	184		.		.

#store gaps of genome assemblies in hash [references]
my $gapgen_ref;
print $logone_fh "\t--- Dealing with...\n";
for (my $f=0;$f<=$#gapfiles;$f++) {
	print $logone_fh "\t    $gapfiles[$f]...\n";
	$gapgen_ref = DelGet::get_gaps($gapgen_ref,$f,$gapfiles[$f]); #access to gaps of a genome = $gapgen_ref->{$f}, with $f=0 for gen1 etc #bug fix 2016 07 27 = send $gapgen_ref to the sub
}
print $logone_fh "\t--- done\n";


###########################################################
# Process genome1 => get names and lengths + nb of random per scaff for genone
###########################################################
print $logone_fh "\n---GET GENOME1 (OUTGROUP) LENGTHS & INFOS----------------------------\n";

#Extract filename
my $genone_name = DelGet::filename($genone);

print $logone_fh "\t--- Getting total length of scaffolds to consider + nb of randomizations per scaffold in hash\n";
close $logone_fh;
#Get total length 
my ($totlength) = DelGet::get_tot_len_filtered($genone,$log,$minlen);
#Getnb of random positions to get per sequence in hash [reference]
my $infofile = "$path/$genone_name.rlen$anch_dist.alen$anch_len.$randtot"."rand.infos.tab";
my ($geninfos_ref,$dbIDs_ref) = DelGet::get_all_len_filtered($genone,$log,$minlen,$totlength,$infofile,$randtot);
open($logone_fh, ">>$log") or die "\t    ERROR - can not open to write log file $log $!\n";
print $logone_fh "\t--- done\n";


######################################################################################################################
# Randomize to get anchors in genome1
######################################################################################################################
#number of rounds [will be re-initialized to last previous round if there is $OKregions file loaded]
my $r = 1;

###########################################################
# DEAL WITH PREVIOUS OUTPUT IF NEEDED
###########################################################
print $logone_fh "\n---PREP FILES---------------------------------------------\n";

# Put previously randomized regions in a file instead, where new ones will be appended too, to avoid issue of buffering
my $prev_reg = "$pathtemp/_OKregions.previous-coords.tab";
close $logone_fh;
unless ($if_OKregions eq "no") {
	$r = DelGet::get_previous_OKreg($path,$r,$OKregions,$log,$c,$prev_reg);
} else {
	open($logone_fh, ">>",$log) or die "\t    ERROR - can not open to write log file $log $!\n";
	print $logone_fh "\t--- NO previous output file has been defined in CONFIG file\n";
	close $logone_fh;
}
open($logone_fh,">>",$log) or die "\t    ERROR - can not open to write log file $log $!\n";


###########################################################
# Prep outputfiles
###########################################################
my $finalout = "$pathtemp/_OKregions.tab";
print $logone_fh "\t--- Prep output files in $pathtemp\n";
close $logone_fh;
open(my $finalout_fh,">",$finalout) or die "\t    ERROR - can not open to write $finalout $!";
print $finalout_fh "#REGION\t$genone\t\t\t\t\t\t\t\t\t$gentwo\t\t\t\t\t\t\t\t\t\t\t\t\t\t$genthree\t\t\t\t\t\t\t\t\t\t\t\t\t\n";
print $finalout_fh "#ID\ttype\tGname\tGstart\tGend\ttype\tGname\tGstart\tGend\tDIST\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\tstrand\tDIST\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\tstrand\tDIST\n\n";
close $finalout_fh ;

my $gen2out = "$pathtemp/_Gen2.dist.tab";
my $gen3out = "$pathtemp/_Gen3.dist.tab";
my ($gen2out_fh,$gen3out_fh);

# output:
#region1     scaffold/chr     start     end     +     Mluc     /mnt/disk3/genomes/Myotis_lucifugus/Myotis_7x.fasta
my $regions_posi = "$pathtemp/_OKregions.posi.tab";
open(my $regions_posi_fh,">",$regions_posi) or die "\t    ERROR - can not open to write $regions_posi $!";
print $regions_posi_fh "#region_name\tchr\tstart\tend\tstrand\tgenome_ID\tgenome_location\n\n";
close $regions_posi_fh;

#intermediate files during loop
my $anchors = "$pathtemp/data_anchors-posi-fa-blatout";
my $blats = "$pathtemp/data_highest-scores";
mkdir ($anchors, 0755) or die "\t    ERROR - Couldn't mkdir $anchors $!";
mkdir ($blats, 0755) or die "\t    ERROR - Couldn't mkdir $blats $!";


###########################################################
# LOOPS
###########################################################
open($logone_fh,">>",$log) or die "\t    ERROR - can not open to write log file $log $!\n";
print $logone_fh "\n---STARTING LOOPS-----------------------------------------\n";
close $logone_fh;
#keep track of what was already pick [file $prev_reg will be loaded in it if previous stuff]
my %alreadyrand = ();
#checking values
my $ok = 0;
my $failed_anchor_ref;
my $failed_checks_ref;
my %inversions = ();
#Now loop
until ($ok == $randnb) { #big loop
	#re open where previous positions are and load them => avoid overlaps [hash ref]
	open($logone_fh,">>",$log) or die "\t    ERROR - can not open to write log file $log $!\n";
	print $logone_fh "\n\t---LOADING PREVIOUS OK REG----------------------------\n" unless (! -e $prev_reg);
	close $logone_fh;
	my $alreadyrand_full_ref = DelGet::load_previous_OKreg($prev_reg) unless (! -e $prev_reg);
	
	#keep in mind when something was already randomized to avoid overlaps => reinitialize every time there is randomization
	my %randomized_ends = ();
	
	#now carry on
	open($logone_fh,">>",$log) or die "\t    ERROR - can not open to write log file $log $!\n";
	print $logone_fh "\n---> round $r\n";
	my $tempposi = "$pathtemp/data_anchors-posi-fa-blatout/$r.posi.gen1.tab";
	my $tempposifh;
	print $logone_fh "\t--- getting $a_per_round random positions for this round...\n";
	my $i = 0;
	$failed_anchor_ref->{$r}=0; #bug fix 2015 05 13
	RDMPOSI: until ($i == $a_per_round) {
		#Randomize		
		#(note, @dbIDs here only contains IDs that went through filter - check sub)
		my @dbIDs = @{$dbIDs_ref};
		my $rdmIDposi = int(rand($#dbIDs)); 
		my $rdmID = $dbIDs[$rdmIDposi];
		my $ratio = $geninfos_ref->{$rdmID}[1]; #corresponds to nb of random positions to do here, weighten on length
		my $size = $geninfos_ref->{$rdmID}[0] - ($minlen-2); #-2 bc +1 before, and bc minlen will be added to all values	
		if (($ratio < 1) || ($size < 1)) { #shouldn't happen since now IDs are only the ones of contigs larger than $minlen, but just in case keep this
			$failed_anchor_ref->{$r}++;
			#print $logone_fh "ratio_or_size<1=>next\n";
			next RDMPOSI;
		}	
		my $three_end = int(rand($size)) + 1 + ($minlen-2); #get in base1

		# get other coordinates and store the end
		my $three_start = $three_end - $anch_len + 1; #len of anchor = 100nt
		my $five_end = $three_start - $anch_dist + 1; #len in between is really 300nt
		my $five_start = $five_end - $anch_len + 1; #len of anchor = 100nt
		($randomized_ends{$rdmID})?($randomized_ends{$rdmID} .= ",$three_end"):($randomized_ends{$rdmID} = $three_end);
				
		#3) avoid overlap with any already picked region
		my $ifnextrand = "no";
		# - first, do not overlap current randomizations
		($ifnextrand) = DelGet::check_overlap_reg($randomized_ends{$rdmID},$anch_len,$anch_dist,$three_end,$five_start) if (exists $randomized_ends{$rdmID});
		if ($ifnextrand eq "yes") {
			$failed_anchor_ref->{$r}++;
			next RDMPOSI;
		}				
		# - second, check regions previously stored
		($ifnextrand) = DelGet::check_overlap_reg($alreadyrand_full_ref->{$rdmID},$anch_len,$anch_dist,$three_end,$five_start) if (exists $alreadyrand_full_ref->{$rdmID});
		if ($ifnextrand eq "yes") {
			$failed_anchor_ref->{$r}++;
			next RDMPOSI;
		}
		# => not overlapping with previous stuff => OK to continue 
		
		#4) Checking assembly gaps; start and end storred this time b/c size of gaps may change		
		($ifnextrand) = DelGet::check_overlap_gap($gapgen_ref->{0}->{$rdmID},$three_end,$five_start) if (exists $gapgen_ref->{0}->{$rdmID});
		if ($ifnextrand eq "yes") {
			$failed_anchor_ref->{$r}++;
			next RDMPOSI;
		}
				
		#everything OK => print in posifile 
		print $logone_fh "\t    => $i regions, OK\n" if ($i == $a_per_round);
		open($tempposifh,">>$tempposi") or die "\t    ERROR - can not open to write $tempposi $!";
		print $tempposifh "reg".$r."-".$i."#5"."\t$rdmID\t$five_start\t$five_end\n";
		print $tempposifh "reg".$r."-".$i."#3"."\t$rdmID\t$three_start\t$three_end\n";
		$i++;
		close $tempposifh;
		#since there was 1 random succesful for this one, remove 1 #Moved here bug fix 2016 07 27
		$geninfos_ref->{$rdmID}->[1]--;
	}
	print $logone_fh "\t\tNB, failed randomizations round $r = $failed_anchor_ref->{$r}\n";
	print $logone_fh "\t--- done\n";
	
	######################################################################################################################
	# Extract sequences of randomized anchors
	######################################################################################################################
	print $logone_fh "\n\t---EXTRACT SEQUENCES----------------------------------\n";
	my $anchors = "$pathtemp/data_anchors-posi-fa-blatout/$r.anchors.gen1.fa";
	close $logone_fh;
	my $gen1Infos_ref = DelGet::extract_sequences($tempposi,$anchors,$log,$genone);;
	open($logone_fh, ">>$log") or die "\t    ERROR - can not create log file $log $!\n";
	print $logone_fh "\t--- done\n";
	

	######################################################################################################################
	# BLAT this fasta file against genome2 and genome3
	######################################################################################################################
	print $logone_fh "\n\t---BLAT-----------------------------------------------\n";
	print $logone_fh "\t    blat in progress, against gen2 = $gentwo\n";
	my $blatout = "$anchors.blat.gen2.out";
	system "$BLATSOFT -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 $gentwo $anchors $blatout";
	print $logone_fh "\t    blat in progress, against gen3 = $genthree\n";
	my $blatout2 = "$anchors.blat.gen3.out";
	system "$BLATSOFT -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 $genthree $anchors $blatout2";
	my @blatfiles = ($blatout,$blatout2);
	print $logone_fh "\t--- done\n";

	
	######################################################################################################################
	# PARSE BLAT OUTPUTS
	######################################################################################################################
	print $logone_fh "\n\t---PARSE BLAT-----------------------------------------\n";
	#get highest hit score for each query + keep info of how higher it is (if relevant)
	my @HSfiles = ();
	tie @HSfiles, 'Array::Unique';
	print $logone_fh "\t--- getting highest scores from blat outputs...\n";
	close $logone_fh;
	for (my $i = 0; $i<=$#blatfiles; $i ++) {
		my @HSfiles_temp = DelGet::parse_blat($log,$blatfiles[$i],$genIDs[$i],$pathtemp);
		push (@HSfiles,@HSfiles_temp);
	}
	open($logone_fh, ">>$log") or die "\t    ERROR - can not create log file $log $!\n";
	print $logone_fh "\t--- done\n";
	#print $logone_fh "\t--- files to process = @HSfiles\n";
	
	######################################################################################################################
	# CHECK OUTPUTS, one per region
	######################################################################################################################
	print $logone_fh "\n\t---CHECKING HIGHEST SCORES----------------------------\n";
	print $logone_fh "\t    [if hit, if same target, if score is high enough, if same strand, if overlap gaps in assembly...]\n";
	#loop on region
	REGIONS: for (my $i = 0; $i<=$#HSfiles; $i ++) {
		my $HSfile = $HSfiles[$i];
		my %gentwo;
		my %genthree;
		my @anchor_names;
		open (HS,"<$HSfile") or die print $logone_fh "\t    ERROR - can not open file $HSfile $!";
		while (<HS>) { 
			chomp (my $line = $_);
			my @line = split(/\s+/,$line);
			my ($strand,$t_name,$t_start,$t_end,$type,$gen,$ratio,$score) = ($line[8],$line[13],$line[15],$line[16],$line[22],$line[23],$line[27],$line[25]);
			push(@anchor_names,$line[9]);
			#			0		1		2		 3		4	 5		6
			my @vals = ($strand,$t_name,$t_start,$t_end,$gen,$ratio,$score);
			if ($gen eq $IDgen2) {
				$gentwo{$type} = \@vals;
			} elsif ($gen eq $IDgen3) {
				$genthree{$type} = \@vals;
			}
		}
		close HS;
		
		print $logone_fh "\t--- parsing $HSfile in progress...\n";
		
		# If one anchor doesn't have a hit
		unless ((exists $gentwo{5}) && (exists $gentwo{3}) && (exists $genthree{5}) && (exists $genthree{3})) {
			#print $logone_fh "\t    one anchor = not hit in gen2 $IDgen2\n" unless ($gentwo{5} && $gentwo{3});
			#print $logone_fh "\t    one anchor = not hit in gen3 $IDgen3\n" unless ($genthree{5} && $genthree{3});
			$failed_checks_ref->{$r}++;
			next REGIONS;
		}
		
		# check same target
		if (($gentwo{5}->[1] ne $gentwo{3}->[1]) || ($genthree{5}->[1] ne $genthree{3}->[1])) {
			#print $logone_fh "\t    not same target in gen2 $IDgen2\n" if ($gentwo{5}->[1] ne $gentwo{3}->[1]);
			#print $logone_fh "\t    not same target in gen3 $IDgen3\n" if ($genthree{5}->[1] ne $genthree{3}->[1]);
			$failed_checks_ref->{$r}++;
			next REGIONS;
		}
		my $t_name_gen2 = $gentwo{5}->[1];
		my $t_name_gen3 = $genthree{5}->[1];
		
		# Check anchor sizes < $a_max_len (150%)
		if (($gentwo{5}->[3] - $gentwo{5}->[2] +1 > $a_max_len) || ($gentwo{3}->[3] - $gentwo{3}->[2] +1 > $a_max_len) || ($genthree{5}->[3] - $genthree{5}->[2] +1 > $a_max_len) || ($genthree{3}->[3] - $genthree{3}->[2] +1 > $a_max_len)) {
			#print $logone_fh "\t    One anchor is > 150nt\n";
			$failed_checks_ref->{$r}++;
			next REGIONS;
		}
		
		#check highest enough
		if (($gentwo{5}->[5] > $maxratio) || ($gentwo{3}->[5] > $maxratio) || ($genthree{5}->[5] > $maxratio) || ($genthree{3}->[5] > $maxratio)) {
			#print $logone_fh "\t    at least one hit did not have scores high enough\n";
			$failed_checks_ref->{$r}++;
			next REGIONS;
		}

		# check strand
		if (($gentwo{5}->[0] ne $gentwo{3}->[0]) || ($genthree{5}->[0] ne $genthree{3}->[0])) {
			$inversions{$r}++;
			next REGIONS;
		}
		
		# check that distance between anchors is not <0
		my $strand_gen2 = $gentwo{5}->[0];
		my ($undef2,$start_gen2,$end_gen2,$dist_gen2) = DelGet::check_anchor_dist($strand_gen2,%gentwo);
		print $logone_fh "\t\tERROR: strand_gen2 undefined? (= $strand_gen2)\n" if ($undef2 == 1);
		next REGIONS if ($dist_gen2 < 0);
		
		my $strand_gen3 = $genthree{5}->[0];
		my ($undef3,$start_gen3,$end_gen3,$dist_gen3) = DelGet::check_anchor_dist($strand_gen3,%genthree);
		print $logone_fh "\t\tERROR: strand_gen3 undefined? (= $strand_gen3)\n" if ($undef3 == 1);
		next REGIONS if ($dist_gen3 < 0);
		
		#check assembly gaps
		if ($gapgen_ref->{1}->{$t_name_gen2}) {
			my ($ifnext2) = DelGet::check_overlap_gap($gapgen_ref->{1}->{$t_name_gen2},$end_gen2,$start_gen2);
			if ($ifnext2 eq "yes") {
				$failed_checks_ref->{$r}++;
				next REGIONS;
			}
		}
		if ($gapgen_ref->{2}->{$t_name_gen3}) {
			my ($ifnext3) = DelGet::check_overlap_gap($gapgen_ref->{2}->{$t_name_gen3},$end_gen3,$start_gen3);
			if ($ifnext3 eq "yes") {
				$failed_checks_ref->{$r}++;
				next REGIONS;
			}
		}	
	
		######################################################################################################################
		# PRINT FINAL OUTPUT
		######################################################################################################################
		print $logone_fh "\t    Went through filters => printing region info in output file\n";
		#				0		1		2		 3		4	 5		6
		#	my @vals = ($strand,$t_name,$t_start,$t_end,$gen,$ratio,$score);
		
		# 1) Regions coords in: Genome1
		my ($gen1_5,$gen1_3,$RID,$regname,$Rtype,$Rstart,$Rend);
		my %regiongen1 = ();
		foreach my $anchor_name (@anchor_names) {
			($RID,$Rstart,$Rend) = split(/\t/,$gen1Infos_ref->{$anchor_name});
			($regname,$Rtype) = split("#",$anchor_name);
			$gen1_5 = "5\t$RID\t$Rstart\t$Rend" if ($Rtype == 5);
			$gen1_3 = "3\t$RID\t$Rstart\t$Rend" if ($Rtype == 3);
			$regiongen1{5} = $Rstart if ($Rtype == 5);
			$regiongen1{3} = $Rend if ($Rtype == 3);
		}
		open($finalout_fh,">>",$finalout) or die "\t    ERROR - can not open to write $finalout $!";
		print $finalout_fh "$regname\t$gen1_5\t$gen1_3\t$anch_dist\t";
		open($regions_posi_fh,">>",$regions_posi) or die "\t    ERROR - can not open to write $regions_posi $!";
		print $regions_posi_fh "$regname\t$RID\t$regiongen1{5}\t$regiongen1{3}\t+\t$IDgen1\t$genone_a\t.\n";
		open($gen2out_fh,">>",$gen2out) or die "\t    ERROR - can not open to write $gen2out $!";
		open($gen3out_fh,">>",$gen3out) or die "\t    ERROR - can not open to write $gen3out $!";
		
		# 2) Regions coords in: Genome2
		my $dist2;
		if ($strand_gen2 eq "+") {
			$dist2 = $gentwo{3}->[2] - $gentwo{5}->[3] + 1;
			print $finalout_fh "5\t$gentwo{5}->[1]\t$gentwo{5}->[2]\t$gentwo{5}->[3]\t$gentwo{5}->[6]\t$gentwo{5}->[5]\t3\t$gentwo{3}->[1]\t$gentwo{3}->[2]\t$gentwo{3}->[3]\t$gentwo{3}->[6]\t$gentwo{3}->[5]\t$strand_gen2\t$dist2\t";
			print $gen2out_fh "$dist2\n";
			print $regions_posi_fh "$regname\t$gentwo{5}->[1]\t$gentwo{5}->[2]\t$gentwo{3}->[3]\t$strand_gen2\t$IDgen2\t$gentwo_a\t.\n";
		} else {
			$dist2 = $gentwo{5}->[2] - $gentwo{3}->[3] + 1;
			print $finalout_fh "3\t$gentwo{3}->[1]\t$gentwo{3}->[2]\t$gentwo{3}->[3]\t$gentwo{3}->[6]\t$gentwo{3}->[5]\t5\t$gentwo{5}->[1]\t$gentwo{5}->[2]\t$gentwo{5}->[3]\t$gentwo{5}->[6]\t$gentwo{5}->[5]\t$strand_gen2\t$dist2\t";
			print $gen2out_fh "$dist2\n";
			print $regions_posi_fh "$regname\t$gentwo{5}->[1]\t$gentwo{3}->[2]\t$gentwo{5}->[3]\t$strand_gen2\t$IDgen2\t$gentwo_a\t.\n";
		}
		# 3) Regions coords in: Genome3
		my $dist3;
		if ($strand_gen3 eq "+") {
			$dist3 = $genthree{3}->[2] - $genthree{5}->[3] + 1;
			print $finalout_fh "5\t$genthree{5}->[1]\t$genthree{5}->[2]\t$genthree{5}->[3]\t$genthree{5}->[6]\t$genthree{5}->[5]\t3\t$genthree{3}->[1]\t$genthree{3}->[2]\t$genthree{3}->[3]\t$genthree{3}->[6]\t$genthree{3}->[5]\t$strand_gen2\t$dist3\n";
			print $gen3out_fh "$dist3\n";
			print $regions_posi_fh "$regname\t$genthree{5}->[1]\t$genthree{5}->[2]\t$genthree{3}->[3]\t$strand_gen3\t$IDgen3\t$genthree_a\t.\n";
		} else {
			$dist3 = $genthree{5}->[2] - $genthree{3}->[3] + 1;
			print $finalout_fh "3\t$genthree{3}->[1]\t$genthree{3}->[2]\t$genthree{3}->[3]\t$genthree{3}->[6]\t$genthree{3}->[5]\t5\t$genthree{5}->[1]\t$genthree{5}->[2]\t$genthree{5}->[3]\t$genthree{5}->[6]\t$genthree{5}->[5]\t$strand_gen2\t$dist3\n";
			print $gen3out_fh "$dist3\n";
			print $regions_posi_fh "$regname\t$genthree{5}->[1]\t$genthree{3}->[2]\t$genthree{5}->[3]\t$strand_gen3\t$IDgen3\t$genthree_a\t.\n";
		}
		$ok ++;
		
		#save that this reg is OK (gen1) to avoid picking any other overlapping one in case there is another loop
		($alreadyrand_full_ref->{$RID})?($alreadyrand_full_ref->{$RID}.=",$regiongen1{3}"):($alreadyrand_full_ref->{$RID}=$regiongen1{3});
		if ($ok == $randnb) {
			print $logone_fh "\t    $ok OK regions [= $randnb => exit]\n";
			print $logone_fh "\t\tNB, failed check after blats round $r = $failed_checks_ref->{$r}\n" if $failed_checks_ref->{$r};
			print $logone_fh "\t\tNB, potential inversions (anchors not same strand) of round $r = $inversions{$r}\n" if $inversions{$r};
			exit;
		}
		close $finalout_fh;
		close $regions_posi_fh;
		close $gen2out_fh;
		close $gen3out_fh;
	} #end loop REGION	
	print $logone_fh "\t    $ok OK regions [still < $randnb => going next round]\n";
	print $logone_fh "\t\tNB, failed check after blats round $r = $failed_checks_ref->{$r}\n" if $failed_checks_ref->{$r};
	print $logone_fh "\t\tNB, potential inversions (anchors not same strand) of round $r = $inversions{$r}\n" if $inversions{$r};
	$r++;
	
	#rewrite file of previous regions, using alreadyrand_full, since this hash is reinitialized and re loaded from file at each run of loop RANDOM
	open(my $prev_out_fh,">$prev_reg") or die "\t    ERROR - can not open to write file $prev_reg $!";
	foreach my $prev_chr (sort keys %{$alreadyrand_full_ref}) {
		my @already = split(",",$alreadyrand_full_ref->{$prev_chr});
		foreach my $prev_end (@already) {
			print $prev_out_fh "$prev_chr\t$prev_end\n";
		}
	}
	close $prev_out_fh;
	close $logone_fh;
} #end loop RANDOM	
######################################################################################################################
close $regions_posi_fh;
close $gen2out_fh;
close $gen3out_fh;
print $logone_fh "\n---DONE---------------------------------------------------
 -> $ok regions
 -> see outputs in $pathtemp/ 
----------------------------------------------------------\n\n";
close $logone_fh;
exit;



#STRUCTURE OF HighestScores.tab
#0	1	2	3	4	5	6	7	8	9		10	11	12	13							14		15		16		17	18		19	20			21		22	23		24		25	26				27
#97	4	0	0	0	0	0	0	+	reg2#3	101	0	101	gi|432108408|gb|KB104502.1|	6866958	2472779	2472880	1	101,	0,	2472779,	reg2	3	gen2	score=	100	ratio_with_2nd=	0.2277



#STRUCTURE OF PSL BLAT OUTPUT AND VARIABLE NAMES TO CALL VALUES
#FROM http://doc.bioperl.org/releases/bioperl-1.4/Bio/SearchIO/psl.html
#0			1			2				3		4				5				6				7				8			9		10			11			12		13		 14			15			16		17	
#$matches  $mismatches  $rep_matches  $n_count  $q_num_insert  $q_base_insert  $t_num_insert   $t_base_insert   $strand   $q_name   $q_length   $q_start  $q_end   $t_name   $t_length  $t_start   $t_end   $block_count  
#18				 19			 20
#$block_sizes    $q_starts   $t_starts

#0		1		2		3	4		5		6		7		8		9			10		11		12	13			14		15		16	17		18			19		 20
# match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	# match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count			
#---------------------------------------------------------------------------------------------------------------------------------------------------------------																	
# 17	0	0	0	0	0	0	0	+	GL429775_6626316_6626416_300bp_seq2_luci5'	101	83	100	gi|432119242|gb|KB098985.1|	2558592	1429419	1429436	1	17,	83,	1429419,



