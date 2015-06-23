#!/usr/bin/perl -w

##########################################################
# Author  :  Aurelie Kapusta
# version :  4.2 (see updates)
# email   :  4urelie.k@gmail.com
# PURPOSE :  General question at the base of writing this pipeline = assess medium size deletions between species.
#			 See lower ($usage and $help) for details.
#			 Written for 3 species but can be adapted to more quite easily - besides gap analysis that will be harder
##########################################################
# REQUIRES :
#	DelGet--0--CONFIG.pl => text file where you can edit paths and some variables
#	DelGet--1--get_regions.pl
#		=> requires UCSC BLAT stand alone software, see http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3
#	DelGet--2--get-sev-seq_align.pl
#       => requires MUSCLE (http://www.drive5.com/muscle/downloads.htm) 
#       => depending on region sizes, may require as well KALIGN (http://msa.sbc.su.se/cgi-bin/msa.cgi?mode=downloads)
#	DelGet--3--gapfreq.pl
##########################################################
use strict;
use warnings;
use Getopt::Long;
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
my $version = "4.3";

my $changelog = "
# UPDATES:
#	- v1.0 = 02 May 2013
#	- v2.0 = 14 May 2013
#		lots of bug fixing beetween v1.0 and v2.0
#		changed to make it a pipeline to launch the 3 scripts, so that only pipeline script is changed when genomes and parameters are changed
#	- v2.1 = 20 May 2013
#		check that distance between anchors >=0
#	- v2.2 = 26 Jun 2013
#		possibility of loading previous output file, to avoid overlaps => append new regions to a previous run
#	- v3.0 = 14 Oct 2013
#		see other scripts for more details about modifications
#		In this new version, runs for anchors are made 100 by 100 only => change on loops and pieces of scripts are moved from --1-- to this one.
#		This is because there were issues of script stopping to get regions, different numbers based on organisms
#		So now that it loops on smaller amounts, integration of the other scripts to the pipeline is possible => get more regions and have analysis done \"automatically\"
#	- v3.2 = 18 Oct 2013
#		added Kalign software to align big regions
#		replaced all \"muscle\" in names by \"align\"	
#	- v3.3 = 05 Nov 2013
#		\$oktot was not incrementing [muscle in file name instead of align]
#		Regions were still overlapping => see Deletions_pipeline--1--get_regions_ak.pl
#	- v3.4 = 08 Nov 2013
#		changes in Deletions_pipeline--3--gapfreq_ak.pl and Deletions_pipeline--3--gapfreq_lc-uc_ak.pl scripts
#			=> skip the A in situations like TGC-----A--TGC, e.g. gap opening is - after the C and gap ending is - before the T. But gap len is still good.
#		change in Deletions_pipeline--3--gapfreq_lc-uc_ak.pl script
#			=> gap is considered \"in masked\" if 75% of the gap is lower case in at least one other species (it was 85% before)
#	- v3.5 = 14 Nov 2013
#		Regions were still overlapping, issue was DURING run [wrong coordinate stored] => see Deletions_pipeline--1--get_regions_ak.pl
#	- v3.6 = 07 Jan 2014
#		Changes in Deletions_pipeline--1--get_regions_ak.pl to (try to) avoid the running with no printing outputs
#	- v4.0 = 05 Feb 2014
#		Changes to use subroutines => no functional changes but different aspect of the code
#	- v4.1 = 04 Mar 2014
#		Little bugs corrected, in
#			Deletions_pipeline--1--get_regions_ak.pl (getting previous regions)
#			Deletions_pipeline--3--gapfreq_cat-outputs_ak.pl (errors in concat, see script)
#			fasta-aln_RW-output-with-lowcases_from-RMout_ak.pl (files named align and not muscle now)
#   - v4.2 = 27-29 Jan 2015
#		BEGIN stuff + added some die checks
#   - v4.3 = 20 Mar 2015
#		chlog and few details like that

#   TO DO
#		- config file printed by the pipeline + read by a subroutine instead
#		- merge everything in one script that uses modules -> will be easier to kill...
#       - also would allow to not randomize and load a set of coordinates => interrogate specifically some regions
#       - make is possible to use mor than 3 species as long as a tree is provided
#       - parallelize when extract and align
#       - print STDERR and not a log
\n";

#################################################################
# Usage and help
#################################################################
my $usage = "\nUSAGE [v$version]:
	perl script.pl <config_file>

	Typically:
		perl DelGet.pl DelGet--0--CONFIG_manipname.pl

	Please type:
	\"perl DelGet.pl -h\" for more details about the pipeline and config file
	\"perl DelGet.pl -chlog\" for more details about the various updates
	\"perl DelGet.pl -version\" for the version number\n\n";
	
my $help = "
PURPOSE :   General question = get medium size deletion rates
            Requires the UCSC BLAT stand alone software, see http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3
            Written for 3 species (referred as 1, 2 and 3)
                1) From a \"reference\" genome (=1), get random regions as: [anchor1]=====Genome1.region=====[anchor2]
                2) Blat anchor1 and anchor2 sequences agains genome2 and genome3
                3) Get length between anchors: [anchor1]=====Genome2.region=====[anchor2] and [anchor1]=====Genome3.region=====[anchor2]

            NOTE THAT:
                - anchors are chosen based on their score + if second highest score is < XX% of the highest score (see \"maxratio\" variable in CONFIG file)
                - anchors need to be on same target sequence in both genomes
                - anchors on genome2 and genome3 need to be < XX nt (\"a_max_len\" variable in CONFIG file);
                - anchor1 and anchor2 hits need to be on same strand
                - and obviously no assembly gaps in any of the regions
				
				
HOW TO RUN THE SCRIPT:
            perl DelGet.pl <config_file>

            Typically:
                perl DelGet.pl DelGet--0--CONFIG_manipname.pl
	
	
IMPORTANT: CONFIG file 
            This is were you NEED to define the 3 genome locations + gap files, by editing path between the quotes.
                These gap files list assembly gaps coordinates. Files need to be as UCSC format, 
                   and obvisouly name of sequences need to match names of sequences in your genome file. 
                If no gap file is available on UCSC for your assembly, use the provided fasta_get_gaps.pl
                If you never see in the log file \"gap in this region in gen1 => next random posi\" or \"there are gaps in this region\", 
                   then there is an issue with these files, please check format etc.
            This is were you NEED to define BLAT software location
            This is were you NEED to define the pipeline parameters
            Because of 2 different indexing happen, you need to define different location for genomes to extract sequences from. Sym links work.
		

OUTPUTS :   0) From this pipeline script
             -> _DelGet.log
                You WANT to check it. Steps are detailed, as well as files names, checking steps etc.
	
            1) FROM SCRIPT DelGet--1--get_regions.pl
            INSIDE FOLDER Deletions.\$c (\$c being a counter to avoid erasing previous results since script will be run several times)
             -> _DelGet--1--get_regions.log
                You WANT to check it. Steps are detailed, as well as files names, checking steps etc.
             -> _OKregions.tab
                Contains all coordinates of anchors in the 3 genomes + distance between anchors [see column titles, pretty explicit]
                Will be used to avoind overlaps so don't change its name.
             -> _Gen2.dist.tab and _Gen3.dist.tab
                Contain only whole regions sizes => usable to plot in R directly
             -> _OKregions.posi.tab
                Contains region positions and infos needed to extract and align regions using script DelGet--2--get-sev-seq_align.pl
            INSIDE FOLDERS Deletions.\$c/data_anchors-posi-fa-blatout and Deletions.\$c/data_highest-scores
             -> + files containing anchor sequences, positions, blat outputs, highest scores
            genome index and total length files will be written at genome location, so make sure it is at a location you have writing permissions.
            file GENOME.XXX.rlenYYY.alenZZZ.infos.tabinfos.tab will be created in \$path [see CONFIG file]
                [with XXX = number of total regions wanted, YYY = length between anchors and ZZZ = anchor length] 
			
            2) FROM SCRIPT DelGet--2--get-sev-seq_align.pl
            INSIDE FOLDER _OKregions.posi.tab.ExtractAlign
             -> region sequences, .fa and .align.fa (aligned)
			 
            3) FROM SCRIPT DelGet--3--gapfreq_ak.pl
            INSIDE FOLDER _OKregions.posi.tab.ExtractAlign.\$c (\$c being a counter to avoid erasing previous results if script is run several times)
             -> _IDgen.gaps.bed [see CONFID file for IDgen]
                Contains all gaps: region_name    start_in_aln    end_in_aln    gapID_in_aln    .    strand    gap_len    aln_length
             -> _IDgen.gaps.specific.bed [see CONFIG file for IDgen]
                Contains specific gaps. See script itself for more details.\n\n";


####################################################################################################################
# Check if help / get argument if not -help
####################################################################################################################
my ($ifh,$chlog,$ifV);
GetOptions ('help' => \$ifh,'h' => \$ifh, 'chlog' => \$chlog, 'version' => \$ifV);
die "$help" if($ifh);
die "\n   Version = $version\n\n" if($ifV);
die "$changelog" if($chlog);


my $config_file = shift @ARGV or die "$usage $!\n" ;

# Initialize needed configuration file variables
our $a_max_len;
our $randtot;
our $randnb;
our $anch_dist; 
our $anch_len;
our $maxratio;
our $mingaplen;
our $multip;
our $BLATSOFT;
our $ALNSOFT;
our $RMSOFT;
our $path;
our $genone;
our $gentwo;
our $genthree;
our $genone_a;
our $gentwo_a;
our $genthree_a;
our $IDgen1;
our $IDgen2;
our $IDgen3;
our $gapone;
our $gaptwo;
our $gapthree;
our @mask;
our @masking_folders;

#now load config file where these variables are set
require "$config_file" or die "\t    ERROR - can't open $config_file $!\n";

#Checkexistence of files
die "\t    ERROR - $genone does not exist?\n" if (! -e $genone);
die "\t    ERROR - $gentwo does not exist?\n" if (! -e $gentwo);
die "\t    ERROR - $genthree does not exist?\n" if (! -e $genthree);
die "\t    ERROR - $genone_a does not exist?\n" if (! -e $genone_a);
die "\t    ERROR - $gentwo_a does not exist?\n" if (! -e $gentwo_a);
die "\t    ERROR - $genthree_a does not exist?\n" if (! -e $genthree_a);
die "\t    ERROR - $gapone does not exist?\n" if (! -e $gapone);
die "\t    ERROR - $gaptwo does not exist?\n" if (! -e $gaptwo);
die "\t    ERROR - $gapthree does not exist?\n" if (! -e $gapthree);
foreach my $mask (@mask) { die "\t    ERROR - $mask does not exist?\n" if (! -e $mask); }
die "\t    ERROR - please check value of numeric variables in configuration file\n" if (($a_max_len !~ /[0-9]+/) || ($randtot  !~ /[0-9]+/) || ($randnb  !~ /[0-9]+/) || ($anch_dist !~ /[0-9]+/) || ($anch_len !~ /[0-9]+/) || ($mingaplen !~ /[0-9]+/));

#Check existence of tools
die "\t    ERROR: $BLATSOFT does not exist?\n" if (! -e $BLATSOFT);
die "\t    ERROR: $ALNSOFT does not exist?\n" if (! -e $ALNSOFT);
die "\t    ERROR: $RMSOFT does not exist?\n" if (($RMSOFT ne "") && (! -e $RMSOFT));

`mkdir $path` unless (-e $path);
my $log = "$path/_DelGet.log";
open(LOG, ">$log") or die "\t    ERROR - can not create log file $log $!\n";
####################################################################################################################
# Now do runs 200 by 200; that way can be killed and prev runs will still be done
####################################################################################################################
print "Pipeline started (v$version)... -> see $log for progression\nCheck here for verbose or errors\n";

print LOG "--------------------------------------------------------------------------------------------------------------------
Pipeline (v$version) started with following parameters (see config file);

PART 1 = get regions
   - total number of randomizations = $randtot 
     [they will be done $randnb by $randnb; this is arbitrarily chosen based on when script was stopping to write outputs for some groups of species]
   - randomizations are done in $genone
      with gap data = $gapone
   - from 1 random position => positions for 2 anchors $anch_len nt long, separated by $anch_dist
   - extract sequences
   - blat => $BLATSOFT
   - anchors are blat against:
      $gentwo
         with gap data = $gaptwo
      $genthree
         with gap data = $gapthree
   - Length of N stretch to be considered as a gap = $mingaplen [typically, assembly gaps are 100]
   - When checking blat hits, anchors need to be:
      1) same scaffold
      2) anchor size < 150% input length => $a_max_len 
      3) same strand
      4) highest score + ration of (second score / score) needs to be > $maxratio
      5) region can't overlap with gaps

PART 2 = extract and align
   - with muscle => $ALNSOFT
   - except if regions are > $multip x $anch_dist nt
   - using genome files in specific folder (otherwise index problem), as follow:
      $genone_a
      $gentwo_a
      $genthree_a  
   
PART 3 = mask, if libraries defined in config file
   - with Repeat Masker => $RMSOFT

PART 4 = analyze gaps

--------------------------------------------------------------------------------------------------------------------
NOW RUNNING...
--------------------------------------------------------------------------------------------------------------------\n";

# PIPELINE LOOP
my $oktot = 0;
until ($oktot == $randtot) { #ie total nb defined in config file
	print LOG "$oktot regions are ok\n";
	#Folder for outputs, for each run; need to be provided to each script
	my ($pathtemp,$c) = DelGet::make_out_dir($path);
	print LOG "    Round $c -> folder Deletions.$c\n";
	print "\n--------------------------------------------------------------------------------------------------------------------\nRound $c -> folder Deletions.$c\n";
	
	# 1) call script DelGet--1--get_regions.pl
	print "\n----------------------------------------------------------\nErrors or verbose from DelGet--1--get_regions.pl...\n";
	print LOG "----------------------------------------------------------\n   Running DelGet--1--get_regions.pl...\n";
	system "perl $BIN/DelGet--1--get_regions.pl $config_file $pathtemp";
	print LOG "   ...done\n----------------------------------------------------------\n\n";
	print "\n";
	
	# 2) call script to extract and align
	my $input = "$pathtemp/_OKregions.posi.tab";
	print "\n----------------------------------------------------------\nErrors or verbose from DelGet--2--get-sev-seq_align.pl...\n";
	print LOG "----------------------------------------------------------\nRunning DelGet--2--get-sev-seq_align.pl...\n";
	system "perl $BIN/DelGet--2--get-sev-seq_align.pl $config_file $pathtemp $input";
	my $align_dir = "$pathtemp/_ExtractAlign";
	print LOG "   -> see files in $align_dir\n   ...done\n----------------------------------------------------------\n\n";
	
	# 3) if relevant, mask outputs + rewrite lowercases + analyse gaps
	my @files = `ls $align_dir`;
	unless ($#mask == 0) {
		print LOG "----------------------------------------------------------\nMasking extracted sequences...\n";
		for (my $i = 0; $i <= $#mask; $i++) {
			my $lib = $mask[$i];
			my $rm_log = "run_rm.$masking_folders[$i].log";
			print "\n----------------------------------------------------------\nErrors or verbose from Repeat Masking and gap analysis in masked sequences...\n --- -> with $lib...\n";
			print LOG "     -> with $lib...\n";
			foreach my $f (@files) {
				if ($f =~ /^reg[0-9]+-[0-9]+\.fa$/) {
					system "perl $RMSOFT -lib $lib -e ncbi -xsmall -nolow -noisy $align_dir/$f > $align_dir/$rm_log";
				}
			}
			print LOG "        ...done\n";
		
			mkdir my $dir_out = "$align_dir/$masking_folders[$i]";
			my @masked = `ls $align_dir`;
			if ("@masked" =~ /.*.fa\.masked.*/) { #rewrite, for each masking [if relevant]		
				print LOG "        Rewriting .fa.align.fa files with masked regions as lowercases...\n";
				print "\n --- Verbose and errors for rewriting .fa.align.fa files with masked regions as lowercases...\n";
				system "perl $BIN/fasta-aln_RW-with-lc_from-RMout.pl $align_dir";
			} else {
				print LOG "        No masked sequences for $masking_folders[$i] masking, skip moving files and gap analysis\n\n";
			}
			if ("@masked" =~ /.*.fa\.out.*/) { #moving files even if no masked sequences [when relevant]
				print LOG "        Moving RM output files (in $masking_folders[$i])...\n";
				print "\n --- Errors for moving RM output files (in $masking_folders[$i])...\n";
				system "mv $align_dir/*.fa.out $dir_out/";
				system "mv $align_dir/*.fa.cat* $dir_out/";
				system "mv $align_dir/*.fa.tbl $dir_out/" if ("@masked" =~ /.*.fa\.tbl.*/);
				system "mv $align_dir/*.fa.masked* $dir_out/" if ("@masked" =~ /.*.fa\.masked.*/);
			}
			if ("@masked" =~ /.*.fa\.masked.*/) { #analyze masked, for each masking 
				print "\n --- Errors for analyzing gaps in masked sequences ($masking_folders[$i])...\n";
				print LOG "        Analyzing gaps...\n";
				system "perl $BIN/DelGet--3--gapfreq_lc-uc.pl $dir_out $IDgen1,$IDgen2,$IDgen3 > $dir_out/gapfreq.lcuc.log";
				print LOG "        ...done\n\n";
			}
		}
		print LOG "...masking(s) done\n----------------------------------------------------------\n\n";
		print "\n";
	}	
	
	#analyze all
	print "\n----------------------------------------------------------\nErrors or verbose from gap analysis in ALL sequences...\n";
	unless (-f "$align_dir/gapfreq.log") {
		print LOG "----------------------------------------------------------\nAnalyzing all gaps...\n";
		system "perl $BIN/DelGet--3--gapfreq.pl $align_dir $IDgen1,$IDgen2,$IDgen3 > $align_dir/gapfreq.log";
		print LOG "...done\n----------------------------------------------------------\n\n";
	} else {
		print LOG "----------------------------------------------------------\nAnalysis of all gaps already done\n   ...done\n   ----------------------------------------------------------\n\n";
	}
	#ident $oktot
	my $nb = `ls $align_dir | grep -c .fa.align.fa`;
	chomp $nb;
	$oktot += $nb;
	print LOG "\n--------------------------------------------------------------------------------------------------------------------\nRound $c done, $nb regions are analyzed => $oktot regions total\n--------------------------------------------------------------------------------------------------------------------\n";
	print "\nRound $c done, $nb regions are analyzed\n--------------------------------------------------------------------------------------------------------------------\n => $oktot regions total\n--------------------------------------------------------------------------------------------------------------------\n";
} #END OF PIPELINE LOOP

# Now cat the outputs
print "\nNow concatenating outputs";
print LOG "\nNow concatenating outputs\n   command line = $BIN/Utilities/DelGet--3--gapfreq_cat-outputs.pl $path\n";
system "perl $BIN/Utilities/DelGet--3--gapfreq_cat-outputs.pl $path";
print "...done\n--------------------------------------------------------------------------------------------------------------------\n";
print LOG "...done\n--------------------------------------------------------------------------------------------------------------------\n";
close LOG;
exit;



















