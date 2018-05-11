#!/usr/bin/perl -w
#------------------------------------------------------------------------------
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# email   :  4urelie.k@gmail.com
# PURPOSE :  Assess micro (1-30) and medium size (>30) deletion rates between species
#			 perl DelGet.pl -h
#------------------------------------------------------------------------------
BEGIN{
   #what to do on: kill -s ALRM <pid> so I can check where it is if it stalls
   $SIG{ALRM}  = sub {print STDERR "SIGALRM received\n"; print STDERR Carp::longmess; print "\n";};
   #what to do on ^C
   $SIG{INT}  = sub {print STDERR "SIGINT received\n"; print STDERR "\n\n".Carp::longmess; exit;};
   #add a folder in INC
   #unshift(@INC, "~/bin/BioPerl-1.6.901");
}
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN/Lib");
}
use Array::Unique;
use Array::Transpose::Ragged qw/transpose_ragged/;
use Acme::Tools qw(minus);
	
#------------------------------------------------------------------------------
#--- CHANGELOG & USAGE ---------------------------------------------------------
#------------------------------------------------------------------------------
my $VERSION = "5.3";
my $SCRIPTNAME = "DelGet.pl";

my $CHANGELOG;
set_chlog();
sub set_chlog {
	$CHANGELOG = "
# UPDATES:
#	- v1.0 = Apr-May 2013
#	- v2.0 = 14 May 2013
#		Lots of bug fixing beetween v1.0 and v2.0
#       Notably: check step to avoid extracting (and therefore aligning) if region is too big
#                check step to avoid aligning if region*muscle.fa file exists
#		Changed to make it a pipeline to launch 3 scripts, so that when genomes and parameters are changed,  
#          only pipeline script is changed
#	- v2.1 = 20 May 2013
#		Check that distance between anchors >=0
#	- v2.2 = 26 Jun 2013
#		Possibility of loading previous output file, to avoid overlaps => append new regions to a previous run
#	- v2.3 = 19 Jul 2013
#		Maximum anchor length = 150% input (and not fixed at 120nt)
#       Change of \$multip variable value (dep. if 10kb or 25kb)
#       Even if no sequence extraction, do alignment
#	- v2.4 = 08 Sep 2013
#		Corrected error at the step \"avoid overlap with any already picked region\" 
#          when positions are randomized that caused regions to overlap
#	- v3.0 = 14-17 Oct 2013
#		Change way of getting previous output
#       Re integration into pipeline (use of config file, changes in paths)
#		Runs for anchors are made X by X only => change on loops and pieces of scripts are moved from called script to the pipeline.
#		This is because there were issues of script stopping to get regions, different numbers based on organisms
#		So now that it loops on smaller amounts, integration of the other scripts to the pipeline is possible 
#          => get more regions and have analysis done \"automatically\"
#       Modification of regions extraction step (dealing with chr names)
#       Muscle parameter change for very long sequences, to avoid error *** ERROR ***  MSA::GetLetter
#	- v3.2 = 18 Oct 2013
#		Added Kalign software to align big regions (muscle still not working)
#		Replaced all \"muscle\" in names by \"align\"	
#	- v3.3 = 05 Nov 2013
#		\$OKTOT was not incrementing, because \"muscle\"  in file name instead of \"align\"	
#		Regions were still overlapping (\%alreadyrand_full was not declared properly)
#       Added version follow up
#	- v3.4 = 08 Nov 2013
#		Changes in gap analysis:
#		   Skip the A in situations like TGC-------A---TGC, e.g. gap opening is - after the C and gap ending is - before the T. 
#             But gap len is still good.
#		   Gap is considered \"in masked\" if 75% of the gap is lower case in at least one other species (it was 85% before)
#	- v3.5 = 14 Nov 2013
#		Regions were still overlapping, issue was DURING run: wrong value was stored (not end of region)
#       + during debugging I ended up changing ways of checking for overlap and double checked that assembly gap overlap was not happening
#	- v3.6 = 07 Jan 2014
#		Try to avoid the running with no printing outputs
#		   => All exists tests -> define tests
#		   => store of previous positions to avoid in a doc that will be close and reopened at each round
#	- v4.0 = 05 Feb 2014
#		Changes to use subroutines => no functional changes but different aspect of the code
#	- v4.1 = 04 Mar 2014
#		Little bugs corrected:
#		   In getting previuous regions (when Deletions.2, was trying to open Deletions.0/_OKregions.tab instead of in Deletions.1)
#		   Errors in concat during gap analysis
#		   Issue with rewrite outputs with lower cases for masked alignments (files named align and not muscle now)
#   - v4.2 = 27-29 Jan 2015
#		for Github first upload: @INC stuff in the BEGIN statement + added some die checks       
#   - v4.3 = 20 Mar 2015
#		Chlog and few details like that
#   - v4.4 = 13 May 2015
#       Few changes of no consequences while debugging
#		In sub get_all_len_filtered, if ratio = 0, put it at 1 - otherwise it can never go through when many sequences and few are asked
#       Set counter \$failed_anchor_ref->{\$r} to 0
#   - v4.5 = 22-25 Jun 2015
#       Kalign command line issue - options were different => update
#		Change in gap analysis, so that IDs of genome can be assembly IDs without this script returning errors
#   - v4.6 = 26-27 July 2016
#       Bug fix in loading gaps (suroutine get_gaps, only last genome was loaded)
#       Bug fix in decrementing the number of randomizations to do in each sequence, it was done before filters,
#          so the number was going down before a succesfull region was picked for a sequence, and could explain the stalling
#       Just in case having \"our\" variables is an issue and also because it's simply better: 
#          now config file prints and is loaded with a sub routine
#       Introduce a kill switch if more than 5 rounds with 0 OKregions 
#   - v5.0 = 11 Nov 2016
#       Lots of rewriting, to allow different application (for Jainy Thomas):
#		   => load a RM .out file with coordinates to use instead of random positions
#             and use flankings of the coordinates (e.g find internal deletions, could even be solo LTRs 
#             or complete element polymorphism), OR load them 2 by 2 and use flankings of 2 features 
#             (e.g. find illegitimate recombination between simple elements, could also be LTR/LTR recombination)
#       If 3 species => same as DelGet.pl, get oriented gaps. If only fewer or more than 3 species, just list them/align if asked. 
#   - v5.1 = 26 Mar 2018
#       Major update
#       -> Mostly, implement gap rates for more than 3 species.       
#          So now a newick tree is directly loaded, and gap analysis is done with gap coordinates dubtractions
#          such as MAFmicrodel. Bedtools intersectBEd: with -f 0.85 for > 30nt and -f 0.80 for 1-30
#       -> Remerge of all code into this one file => clean up and straighten the code
#   - v5.2 = 04 Apr 2018
#       Bug fix in randomization stallings
#       Bug fix in large files loaded, still need to do the blats with rounds => randnb by randnb
#       Bug fix TOTLEN
#       Use a 1 liner for the fasta 1 line rewriting
#       Fixed formatting of _OKregions.tab [that bug had pretty drastic consequences]
#       Bug fix loading previous _OKregions.cat.tab
#   - v5.3 = 11 May 2018
#       Couple of changes in the kill switch loops [but scripts are still stalling]
#       Re-implemented the change of v 3.4 => Skip the nucleotide in situations like: TGC-------A---TGC,
#          e.g. gap opening is - after the C and gap ending is - before the T (but length is still good)
#       Added a --restart flag, for the code to check the last folder, delete if incomplete before running
#          and print the info about it in the log
\n";
#   TO DO LIST / IDEAS
#       - clean up sub after each round, to delete blat outputs and highscores that didn't pass filters, unless --debug set
#         also add many more err messages if --debug is set
#       - reintegrate the masking possibility, but fix the length of aln: aln with masked regions, not complete total length...!
#       - detect new insertions... write them out if not alrady masked
#       - parallelize when aligning regions (it is quite a limiting step)
#       - load the My from a newick tree as well and apply that to get the rates
#       - call R and make some graphs + write the code...?
}

my $USAGE = "\n USAGE [v$VERSION]:
   Typically:
     nohup perl $SCRIPTNAME DelGet.manipname.txt > DelGet.manipname.log &

   Please type:
     \"perl $SCRIPTNAME -c <config_file_name>\" to print a sample config file that you can edit
	
   Other options:
     -h (BOOL) prints more details about the pipeline and config file
     -l (BOOL) prints change log / updates
     -v (BOOL) prints the pipeline version
\n";

my $HELP;
set_help();
sub set_help {
	$HELP = "
 PURPOSE: 
   Find deletions (using gaps) between orthologous regions or several species

   Regions analyzed can be:
     1) loaded from a file (only Repeat Masker .out files implemented) 
     2) chosen randomly in one of the species (see config files for numbers, length, etc)
   (See below for the steps)
   
 REQUIREMENTS:
   SOFTWARE:
      Blat = http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3
      Muscle = http://www.drive5.com/muscle/downloads.htm
      Kalign = http://msa.sbc.su.se/cgi-bin/msa.cgi?mode=downloads (if large regions to align, muscle crashes)
      BEDtools = https://github.com/arq5x/bedtools2 (Quinlan AR and Hall IM, 2010. Bioinformatics)
   PERL MODULES:
      (already in Lib) Perl Array::Unique = https://metacpan.org/source/SZABGAB/Array-Unique-0.08/lib/Array/Unique.pm
      (already in Lib) Perl Array::Transpose::Ragged = https://metacpan.org/pod/Array::Transpose::Ragged
      Perl Acme::Tools
	
 CITATION:
   Kapusta, Suh and Feschotte (2017) Dynamics of genome size evolution in birds and mammals
   doi: 10.1073/pnas.1616702114 (http://www.pnas.org/content/114/8/E1460.full)
   (version 4.6 used in this paper)
   
   Please see at the end of the help how gaps are analysed; 
   there is a difference between v4.6 and v5+ that will affect deletion rates 
       
 USAGE: 
   All options and paths are set in a config file, that can be obtained by running: 
      perl $SCRIPTNAME --conf configfile.txt
    
   This is were you NEED to define parameters, file paths, etc, by editing (if needed) 
   the lines that do not start with a # and have a = sign
     -> The gap files list assembly gaps coordinates; you may obtain these files from UCSC, 
        let the script make them, or use the script fasta_get_gaps.pl (in Utils).
        If you never see in the log file \"gap in this region in genX => next random posi\" or 
        \"there are gaps in this region\", amd you expect gaps, there may be an issue with these files.

   Once the config file is edited, typically:
     nohup perl DelGet.pl DelGet.manipname.txt &

   Sometimes the code stalls - if no 'recent' Del/Deletion.X folders, 
     and if the only job still running is $SCRIPTNAME DelGet.manipname.txt
     (or if it has been running for too long)
     consider killing the talled job and restarting, with:
     nohup perl $SCRIPTNAME --restart DelGet.manipname.txt &
     The restart flag will check the last folder and delete it if needed.
		
 OUTPUTS: 
    0. Main log file: <path>_DelGet.log
      Steps are detailed, as well as files names, etc.

    1.'Get regions' step 
      (log Delget_1_get-regions.log - check for any errors)
      The genome index and total length files will be written at genome location, 
         so make sure it is at a location you have writing permissions.
      If regions are obtained randomly, a file GENOME.tot-XXX.minYYY.infos.txt will be created in the output directory 
      (path variable in the config file), with XXX = number of total randomizations estimated, YYY = minlen from the config file  
      FOLDER Deletions.X (X being a counter for subfolders)
         -> _OKregions.tab
            All coordinates of anchors in the genomes + distance between anchors
            Will be used to avoid overlaps so do not change its name if you plan to add runs
         -> _GenX.dist.tab and _GenX.dist.tab
            Whole regions sizes => usable to plot in R directly
         -> _OKregions.posi.tab (if the still_align variable is set to yes)
            Region positions and infos needed to extract and align regions
      FOLDERS Deletions.X/data_anchors-posi-fa-blatout and Deletions.X/data_highest-scores
         -> files containing anchor sequences, positions, blat outputs, highest scores
		
    2.'Extract and align' step 
      (log Delget_2_extract-align.log - check for any errors)
      FOLDER _ExtractAlign
         -> region sequences, .fa and .align.fa (aligned)
		 
    3.'Gap analysis' step 
      (log Delget_3_gap-analysis.log - check for any errors)
      FOLDER _GapAnalysis
         -> All gaps in bed files; split by species, 1-30 and >30 nt
      FOLDER Deletions.X
         There will be two outputs (for 1-30 nt and >30 nt) that will contain these columns:
            count	amount	aln_len	%_of_aln	del_nt/10kb
         A rate can be obtained for the branch by dividing 'del_nt/10kb' by the corresponding My.
    
    4.'Gather final data' step
      (log <path>_DelGet.log - check for any errors)
      The steps 1 to 3 will loop until enough alignment is obtained (variable randtot in config file). 
      Then when the pipeline ends, outputs from all Deletions.X folders will be gathers
      and named after the path variable of the config file:
         <path>_1-30_nt.tab
         <path>_30-x_nt.tab
      with similar columns as the intermediate outputs in Deletions.X folders.
      It is also possible to obtain preliminary data using the Util script DelGet_cat-outputs.pl
            
 STEPS: 
    I. Get the regions where deletions are checked: can be random or loaded from a file
       1) If they are random, then a random position is selected in the first species listed
          (if 3 species are specified and if the goal is deletion rates in the most recently diverged species,
          then the outgroup will be the species listed first - see config file) 
          Then a region is defined, using the anch_dist from the config file, such as:
          [anchor1.5']<=====anch_dist=====>[anchor1.3']
       2) If they are loaded from a file, they can be loaded such as:
          - use flankings of the coordinates (e.g to find internal deletions, 
            could even be solo LTRs / complete element polymorphism), 
          - use flankings of 2 consecutive features (e.g. find illegitimate recombination between simple elements)
            sliding window of 1 (any HR event between copies) or 2 (for LTR/LTR recombination, 
            if only LTRs of proviruses with 2 LTRs left are in the input file)
            
    II. Get orthologous regions
       1) the anchors of the regions as defined above (anchor length is defined by anch_len in the config file)
       2) Blat anchor1.5' and anchor1.3' sequences agains the other 2 genomes                       
          NOTE THAT the blat hits will be considered non valid (filtered out) if:
          - the second highest blat score is too close to the highest score (maxratio in the config file)
          - 5' and 3' anchors hits are not same target sequence (same scaffold/contig)
          - blat hit (new anchor length) is > a_max_len nt (defined in the config file) 
          - 5' and 3' anchors hits are not on the same strand
          - if assembly gaps overlap any of the regions
          NOTE that this step may take a while for randomized regions if species are very diverged or if the assemblies 
          are of poor quality, and that if 5 consecutive rounds with X numbers of trials (a_per_round in config file) 
          found no valid orthologous regions the pipeline will exit.

    III. Process. The target regions are as follow: [anchor.5']<=====Genome.region=====>[anchor.3']
       1) get length between the 3' of the 5' anchor and the 5' of the 3' anchor in all the genomes 
          -> _OKregions.tab files in the Deletion.X directories
       2) list them in an output table [concatenation of the OK regions]
       3) if the still_align option is activated: 
          - extract the sequences (with anchors) in the 3 genomes
          - align them
       4) If analyze_gaps option is activated:
          - process gaps
          - concatenate the outputs with the Util script DelGet_cat-outputs.pl
             This script can be called anytime with the directory Del as argument,
             to get preliminary results even if the runs are not done
             
         Notes on gap analysis:
         1) In all available version, unique nucleotides in the middle are skipped, for example:
           	    TGC-------A---TGC
             Here the A is ignored, and the gap opening is - after the C and gap ending is - before the T
             (the length of the total gap will be maintained)
           
         2) Overlapping gaps. When the overlap is < 80% for 1-30 and < 85% for >30nt, they are considered specific to that species.
            In a case like this one:
           	    outgroup ACCGTGTGTATGTGTGTGTGCGTGCGCGCGTATGTGTCTGTCTGTGTGCGTGTCTGTACGTGTATATAT
           	    species1 ACTGTGTGTGTGTGTGTGT------------------------------GTGTCTGTGCGTGTATATAT
           	    species2 ACTGTGTGTATGTGTGTGTG------------TGTGTGTGTGTGTCTGTGTGTCTGTGCGTGTATATAT
           	 With v4.6, species2 has no species-specific gap and species1 has two (1nt and 17nt); the shared portion is ignored (considered as 'shared').
           	 With v5+, both species 'keep' the full length gaps as specific, which is more 'biological'.
           	 If this is a real biological gap, v4.6 would underestimate deletion rates; after comparing same trios of species with v4.6 and v5.3, 
           	 it looks like this mostly affects microdeletion rates.
           	 This could be a technical issue (sequencing of a simple repeat), but since these regions are also more prone to indels, these could be real.
           	 Thus, it is good to keep this type of examples in mind while discussing microdeletion rates.  
\n\n";
}

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

#------------------------------------------------------------------------------
# --- OPTIONS ------------------------------------------------------------------
#------------------------------------------------------------------------------
my ($CONF,$RESTART,$DEBUG,$IFH,$CHLOG,$IFV);
GetOptions ('conf=s'  => \$CONF,
            'restart' => \$RESTART,
            'debug'   => \$DEBUG, 
            'help'    => \$IFH,
            'log'     => \$CHLOG, 
            'version' => \$IFV);

die "$HELP" if ($IFH);
die "\n   $SCRIPTNAME version $VERSION\n\n" if ($IFV);
die "$CHANGELOG" if ($CHLOG);

#Print config file if relevant
print_config() if ($CONF);
exit 1 if ($CONF);

#------------------------------------------------------------------------------
#--- CONFIG FILE & LOAD ALL STUFF ----------------------------------------------
#------------------------------------------------------------------------------
my $CONFIG_FILE = shift @ARGV or die "$USAGE $!\n" ;
die "\t    ERROR - config file $CONFIG_FILE does not exist?\n\n" if (! -e $CONFIG_FILE);

#Load config file
my %CF = ();
load_config();

#Open the main log file & assign STDOUT and STDERR to it (otherwise hard to know where errors are)
my $LOG = $CF{'path'} unless ($CF{'path'} eq ".");
$LOG = $LOG."_DelGet.log";
my $MFH;
#append if $CF{'path'} exists
if (-d $CF{'path'}) {
	open($MFH, ">>", $LOG) or confess "\t    ERROR - can not open to write log file $LOG $!\n"; 
} else {
	open($MFH, ">", $LOG) or confess "\t    ERROR - can not open to write log file $LOG $!\n"; 
}
#local *STDERR = $MFH;

#Check the last folder and delete it if needed
my $IFDELETED;
clean_up() if ($RESTART);
sub clean_up {
	my @f = `ls $CF{'path'} | grep Deletions`;
	foreach my $f (@f) {
		chomp($f);
		if (! -f $CF{'path'}."/".$f."/__1-30.RESULTS.tab") {
			system "rm -Rf $CF{'path'}/$f";
			if ($IFDELETED) {
				$IFDELETED.=",".$f;
			} else {
				$IFDELETED=$f;
			}
		}
	}
	return 1;
}

#Make output directory
if (! -e $CF{'path'}) {
	mkdir ($CF{'path'}, 0755) or confess "\t    ERROR: can not mkdir $CF{'path'} $!";
}

#Check stuff
check_config();
sub check_config {
	#tree file
 	die "\t    ERROR - tree file needs to be defined in the config file\n" if (! $CF{'tree'});
	die "\t    ERROR - tree does not look like the excecpted format, and file $CF{'tree'} does not exist?\n" if ($CF{'tree'} !~ /^\(.*\);$/ && ! -f $CF{'tree'});
	#data folder
	die "\t    ERROR - data folder $CF{'data'} needs to be defined in the config file\n" if (! $CF{'data'});
	die "\t    ERROR - data folder $CF{'data'} does not exist?\n" if (! -d $CF{'data'});	
	#file to load regions if relevant
	if ($CF{'file'}) {
		die "\t    ERROR - $CF{'file'} does not exist?\n" if (! -f $CF{'file'});
		die "\t    ERROR - $CF{'file'} can only be a Repeat Masker output .out file\n" if ($CF{'file'} !~ /\.out$/);
	} else {
	#other variables
		die "\t    ERROR - please check values of numeric variables in configuration file\n" 
			   if (($CF{'a_max_len'} !~ /[0-9]+/) 
				|| ($CF{'randtot'}  !~ /[0-9]+/) 
				|| ($CF{'randnb'}  !~ /[0-9]+/) 
				|| ($CF{'anch_dist'} !~ /[0-9]+/) 
				|| ($CF{'anch_len'} !~ /[0-9]+/) 
				|| ($CF{'mingaplen'} !~ /[0-9]+/));
	}
	
	#Check existence of tools
	die "\t    ERROR: $CF{'blat'} does not exist?\n" if (! -f $CF{'blat'});
	die "\t    ERROR: $CF{'bedtools'} does not exist?\n" if (! -d $CF{'bedtools'});
	die "\t    ERROR: $CF{'alnsoft'} does not exist?\n" if ($CF{'still_align'} && ! -f $CF{'alnsoft'});
	return 1;
}

#print some log, even if appended, just in cases parameters changed
print_log_1();
sub print_log_1 {
	if ($DEBUG) {
		print $MFH "\n-------------------------------------------------------------------------------\n";
		print $MFH "DEBUGGINF MODE:\n";
		print $MFH "   - additional error messages\n";
		print $MFH "   - temp files not deleted\n";
	}
	if ($RESTART && $IFDELETED) {
		print $MFH "\n-------------------------------------------------------------------------------\n";
		print $MFH "RESTARTING PIPELINE:\n";
		print $MFH "   - previous folders were checked,\n"; 
		print $MFH "     and deleted if the __*.RESULTS.tab files were missing\n";
		print $MFH "     Deleted folders = $IFDELETED\n";
	}
	print $MFH "\n-------------------------------------------------------------------------------\n";
	print $MFH "Pipeline (v$VERSION) started with following parameters loaded from $CONFIG_FILE\n";
	print $MFH "I... get regions\n";
	if ($CF{'file'}) {
		print $MFH "     - from $CF{'file'}\n";
	} else {
		print $MFH "     - through randomization:\n";
		print $MFH "        1. get a random position in the outgroup\n";
		print $MFH "        2. get 2 anchors $CF{'anch_len'}nt long, separated by $CF{'anch_dist'} nt\n";
		print $MFH "        They will be done $CF{'randnb'} by $CF{'randnb'} and stop when reach a final alignment length of $CF{'randtot'}nt\n";
	}
	print $MFH "     - extract anchor sequences\n";
	print $MFH "     - blat them against target genomes using $CF{'blat'}\n";
	print $MFH "     - check validity of the hits: anchor hits will be consider non valid if:\n";
	print $MFH "        1. 5' and 3' anchors hits are not same target sequence (same scaffold/contig)\n";
	print $MFH "        2. blat hit (new anchor length) is > 150% input length => $CF{'a_max_len'} nt\n";
	print $MFH "        3. 5' and 3' anchors hits are not on the same strand\n";
	print $MFH "        4. second score / highest score is < $CF{'maxratio'} (the second highest blat score is too close to the highest score)\n";
	print $MFH "        5. if assembly gaps (= gaps of at least $CF{'mingaplen'} nt) overlap any of the regions\n";
	print $MFH "II.. print regions in _OKregions.tab file(s)\n";
	if ($CF{'still_align'}) {
		print $MFH "III. Extract and align\n";
		print $MFH "     - with sofware = $CF{'alnsoft'}\n";
		if ($CF{'alnsoft'} eq "muscle") {
			print $MFH "       Note that Muscle rewrites sequences in upercase; if you see lowercase,\n";
			print $MFH "       it means the alignment was made with Kalign instead (region too large)\n";
		} else {
			print $MFH "       Note that Kalign retains the original case\n";
		}
		print $MFH "     - regions won't be aligned if longer than $CF{'r_max_len'} nt";
		if ($CF{'file'}) {
			print $MFH "\n";	
 		} else {
 			print $MFH " ($CF{'multip'} x $CF{'anch_dist'})\n";	
 		}
#		if ($CF{'mask'}) {
# 			print $MFH "\n IV = mask, if libraries defined in config file\n";
# 			print $MFH "   - with Repeat Masker => $CF{'RMSOFT'}\n";
# 		}
		if ($CF{'analyze_gaps'}) {
			print $MFH "IV.. analyze gaps\n";
			print $MFH "     - Extract gap coordinates\n";
			print $MFH "     - intersect/subtract to determine shared vs specific gaps\n";
		}
	}
	print $MFH "\n";
	print $MFH "-------------------------------------------------------------------------------\n";
	print $MFH "PREP STEPS:\n";
	return 1;
} 

#------------------------------------------------------------------------------
#--- LOAD & PARSE TREE ---------------------------------------------------------
#------------------------------------------------------------------------------
#LOAD & PARSE TREE 
load_newick_tree();

my %NWCK = ();
my $ITER;
my $LEAF = 0;
my %SPIDS = ();
my $MAX = 0;
until ($ITER) {
	parse_newick_tree();
}
my %GROUP = ();
get_groups_from_ids();

my %PAR = ();
get_group_parents();

#Species to consider:
my @SPIDS = @{$SPIDS{$MAX}}; #final list of spIDs
my $NBSP = scalar(@SPIDS); #number of species
print $MFH "     => $NBSP species = @SPIDS\n";
my $OG;
if ($NBSP < 3) {
	$OG = $CF{'file'};
	$OG =~ s/^(.+?)\.out/$1/;
} else {
	$OG = $NWCK{$MAX}; #will have the outgroup + the last group
	$OG =~ s/,[0-9]+//;
	$OG =~ s/[0-9]+,//;
}	
print $MFH "        With outgroup = $OG\n";
die "\t    ERROR - $OG is not a unique spcecies...?\n\n" if ($OG =~ /(\.|\,|\(|\))/);

#------------------------------------------------------------------------------
#--- GENOME AND GAP FILES: LISTS AND STUFF -------------------------------------
#------------------------------------------------------------------------------
print $MFH "\n --- Dealing with the genome files...\n";
my %DATA = (); 
#NB:
#	$DATA{$sp}{'f'} = path/sp.fa
#	$DATA{$sp}{'n'} = sp.fa
#	$DATA{$sp}{'e'} = path/sp.to-extract.fa #just a ln -s of sp.fa
#	$DATA{$sp}{'g'} = path/sp.fa.gaps.min50.tab (typically)
get_genome_files_list();
print $MFH "     ...done\n";

my %INFOS = (); 
my @DBIDS = (); #list of sequences from the outgroup
unless ($CF{'file'}) {
	print $MFH " --- Dealing with the outgroup species $OG\n";
	print $MFH "     -> get sequence lengths and estimate number of randomizations per scaffold\n";
	#NB:
	#	$DATA{$MAX}{'l'} = total length of seqs > the minlength in nt of the outgroup
	#	$INFOS{$id} = ($len,$ratio) #for the outgroup
	get_fa_infos();
	print $MFH "     ...done\n";	
}

print $MFH " --- Dealing with gap files...\n";
get_gap_files_list();
#Load the gap coordinates
my %GAPS = ();
load_gaps();
print $MFH "     ...done\n";	

my %COORDS = ();
if ($CF{'file'}) {
	print $MFH " --- Loading coordinates from $CF{'file'}\n";
	load_file();
	print $MFH "     ...done\n";	
}

#------------------------------------------------------------------------------
#--- MAIN: PIPELINE LOOP -------------------------------------------------------
#------------------------------------------------------------------------------
print $MFH "\n-------------------------------------------------------------------------------\n";
print $MFH "NOW RUNNING:\n";

#MAIN LOOP VARIABLES:
my $KILLSWITCH = 0; #Fail safe
my $OKTOT = 0;      #total amount of nt in alignments, or number, that are done
my $CURRDEL;        #Deletion.X directory
my $C = 0;          #counter of loops

#VARIABLE FOR 'GET REGIONS':
my $R; #counter for region IDs for a given run
my ($ANCHORS,$BLATS); #directories
my ($OKREG,$PREVREG,$POSI,$RLOG,$RFH); #files
my %ALREADYRAND = (); #keep track of what was already picked
my $OK; #total OK regions
my $FAILED_ANCHOR = 0; #reset for each iterations
my %OGINFOS = (); #outgroup infos
my @HSFILES = (); #blat highest scores files
my $FAILED_CHECKS = 0; #reset for each iterations

#VARIABLE FOR 'EXTRACT ALIGN':
my $ALIGN; #directory
my ($ALOG,$AFH); #files
my %KALIGN = (); #to align some regions with kalign if too long

#VARIABLE FOR 'ANALYZE GAPSS':
my $GAPS; #directory
my ($GLOG,$GFH); #files
my %TOTLEN = (); #total lenghth of aln processed per species (per round)
my @GAPFILES = (); #list that will contain the files to process (per round)
my ($ANC,$CURR); #for file names

#MAIN LOOP: will run as long as not enough nt in alignments or if killswitch
until ($OKTOT >= $CF{'randtot'} || $KILLSWITCH == 5) { 
	main_loop();
	if ($KILLSWITCH == 5) {
		print $MFH " => EXITING even if total of ok regions ($OKTOT) < number set in config file ($CF{'randtot'})\n";
		print $MFH "    because $KILLSWITCH rounds without any OK regions.\n";
		print $MFH "    Check for errors, or maybe species are too diverged?\n";
	}	
}
sub main_loop {
	local *STDERR = $MFH;
	print $MFH "   $OKTOT nt in total length of alignments so far\n" if (! $CF{'file'} && $C > 0);	
	my $nb = 0;
	
	#Folder for outputs, for each run
	make_out_dirs(); #where $C is set
	print $MFH "-------------------------------------------------------------------------------\n" unless ($C == 0);	
	print $MFH "MAIN ROUND $C -> folder Deletions.$C ($CURRDEL)\n";

	#1) Get all regions => print all stuff in $CURRDEL
	print $MFH " --- Obtaining regions\n";
	
	#prep output files (print headers)
	$OKREG   = "$CURRDEL/_OKregions.tab";
	$PREVREG = "$CURRDEL/_OKregions.previous-coords.tab";
	$POSI    = "$CURRDEL/_OKregions.posi.tab" if ($CF{'still_align'});
	print $MFH " --- Prep output files\n";
	prep_curr_out_files();
	
	#Put previously randomized regions in a file, where new ones will be appended too
	$R = 1; 
	if ($CF{'if_OKregions'} eq "no") {
		if ($CF{'file'}) {	
			print $MFH "       NOTE: No previous regions will be loaded (expected if anchors loaded from file)\n";
		} else {
			print $MFH "       WARN: No previous regions will EVER be loaded (regions won't necessarily be independent)\n";
		}
	} else {
		if ($CF{'OKregions'} =~ /\w/) {
			print $MFH "       Previous outputs are set to be loaded (file = $CF{'OKregions'})\n";
			get_previous_OKreg(); #will get the last $R if relevant
		}	
	}

	#Get regions
	$OK = 0;
	%OGINFOS = ();
	$RLOG = "$CURRDEL/Delget_1_get-regions.log";
	print $MFH " --- Now looping to get regions => check log in:\n";
	print $MFH "     $RLOG\n";
	open($RFH,">>",$RLOG) or confess "\t    ERROR - can not open to write log file $RLOG $!\n";
	local *STDERR = $RFH;
	my $rkswitch = 0;
	until ($OK >= $CF{'randnb'}) { #safety loop exit if 3 rounds with no anchor that will also ++ the main killswitch
		print $RFH "-------------------------------------------------------------------------------\n";
		print $RFH " ROUND $C - iteration $R\n";
		print $RFH "-------------------------------------------------------------------------------\n";
		#anchors from the outgroup:
		$DATA{$OG}{'posi'} = "$ANCHORS/$R.posi.$OG.tab";
		$DATA{$OG}{'anchors'} = "$ANCHORS/$R.anchors.$OG.fa";	
		$FAILED_ANCHOR = 0;
		if ($CF{'file'}) {
			#read from %COORDS and print in file
			get_regions_from_file();
			$OK=$CF{'randnb'};
		} else {
			#Load previous regions:
			load_previous_OKreg($PREVREG) if (-e $PREVREG); #load them in %ALREADYRAND
			#get new regions
			get_regions();
		}
		if ($CF{'file'} && ! -f $POSI) {
			print $MFH "     WARN: No region went through - exiting\n";
			return 1;
		}	
		extract_anchor_seq();
		blat_anchors();
		@HSFILES = ();
		tie @HSFILES, 'Array::Unique';
		parse_blat();
		$FAILED_CHECKS = 0;
		check_and_print(); #this is where $OK is incremented
		#Now check if should force exit
		if (! $CF{'file'}) {
			$rkswitch++ if ($FAILED_CHECKS == $CF{'a_per_round'});
			if ($rkswitch == 3) {
				print $RFH " EXITING iteration $R (round $C) of randomizing regions, because none could be obtained in the last 3 iterations\n";
				$KILLSWITCH++;
				$R++;
				return 1;
			} else {
				$R++;
				next; #no need to do next steps
			}
		}		
	}	
	close $RFH;
	
	#Extract and align
	if ($CF{'still_align'}) {
		$ALOG = "$CURRDEL/Delget_2_extract-align.log";
		print $MFH " --- Now extracting and aligning regions => check log in:\n";
		print $MFH "     $ALOG\n";
		open($AFH,">",$ALOG) or confess "\t    ERROR - can not open to write log file $ALOG $!\n";
		local *STDERR = $AFH;
		my ($reg,$noaln) = extract_regions();
		my $cnoaln = keys %{$noaln};
		print $MFH "        WARN: $cnoaln regions won't be aligned (check $ALOG)\n" unless ($cnoaln == 0);
		align_regions($reg,$noaln);
		$nb = `ls $ALIGN | grep -c .fa.align.fa` 
	}
	close $AFH;
	
	#Analyze gaps
	if ($CF{'still_align'} && $CF{'analyze_gaps'}) {
		%TOTLEN = ();
		$GLOG = "$CURRDEL/Delget_3_analyze-gaps.log";
		print $MFH " --- Now analyzing gaps => check log in:\n";
		print $MFH "     $GLOG\n";
		open($GFH,">>",$GLOG) or confess "\t    ERROR - can not open to write log file $RLOG $!\n";
		local *STDERR = $GFH;
		@GAPFILES = ();
		analyze_gaps();
		$nb=$TOTLEN{$OG};
		close $GFH;
	}	
	
	#Increment counter and print log	
	chomp($nb);
	print $MFH "ROUND $C DONE\n";
	if ($CF{'file'}) {
		$nb=$CF{'randnb'};
		my $failed = $FAILED_ANCHOR+$FAILED_CHECKS;
		print $MFH "   => $nb regions loaded and checked ($failed did not go through, see $RLOG)\n"; 
	} else {
		print $MFH "   => $nb regions aligned\n" if ($CF{'still_align'} && ! $CF{'analyze_gaps'});
		print $MFH "   => $nb nt of alignment analyzed this round\n" if ($CF{'still_align'} && $CF{'analyze_gaps'});
	}	
	if ($nb == 0) {
		$KILLSWITCH++;
	} else {
		$KILLSWITCH=0; #reset if it went through
	}
	$OKTOT+=$nb;
	return 1;
}

# Now cat the outputs
local *STDERR = $MFH;
if ($CF{'still_align'} && $CF{'analyze_gaps'}) {
	print $MFH "\n-------------------------------------------------------------------------------\n";
	if ($CF{'file'}) {
		print $MFH "DONE - All $CF{'randtot'} features from $CF{'file'} were processed\n"
	} else {
		if ($KILLSWITCH == 5) {
			print $MFH "EXITED - $OKTOT nt from randomization.\n"
		} else {
			print $MFH "DONE - goal ($CF{'randtot'}) reached: $OKTOT nt from randomization.\n"
		}	
	}
	unless ($OKTOT == 0) {
		print $MFH "Concatenate all outputs, with:\n";
		print $MFH "      $BIN/Utils/DelGet_cat-outputs.pl $CF{'path'}\n";
		system "perl $BIN/Utils/DelGet_cat-outputs.pl $CF{'path'}";
	}
}
print $MFH "PIPELINE DONE\n";
print $MFH "-------------------------------------------------------------------------------\n\n";
close $MFH;
exit;


#------------------------------------------------------------------------------
# --- SURBROUTINES -------------------------------------------------------------
#------------------------------------------------------------------------------
sub print_config {	
	die "\n\t    ERROR - $CONF exists, did you mean to run the pipeline with it? If yes, remove -c!\n\n" if (-f $CONF);
	open(my $fh, ">", $CONF) or confess "\t    ERROR: SUB print_config: could not open to write $CONF!\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# CONFIGURATION FILE FOR: DelGet.pl\n";
	print $fh "# VERSION of the pipeline = $VERSION\n";
	print $fh "# if you never used this script, type: perl DelGet.pl -h \n";
	print $fh "#----------------------------------------------------------\n";
	print $fh "# Author  :  Aurelie Kapusta\n";
	print $fh "# email   :  4urelie.k\@gmail.com\n";
	print $fh "# GitHub  :  https://github.com/4ureliek\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# Use this file to edit paths or values anytime you see a line not starting with # and with a =\n";
	print $fh "# do NOT remove the # signs and do NOT put / at the end of the paths; spaces don't matter\n";
	print $fh "\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# DATA VARIABLES - BEHAVIOR OF THE PIPELINE\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# Newick tree:\n";
	print $fh "#   It needs to contain ONLY the species IDs and have at least 3 species to analyze gaps\n";
	print $fh "#   You may provide a file containing the tree, or give the tree on the command line:\n";
	print $fh "    tree = (rheMac3,(hg38,panTro4)); #newick tree, with the ; at the end\n";
	print $fh "#   tree = tree.nwk #file with the tree\n";
	print $fh "\n";	
	print $fh "# Path to genomes & gap files (the files can be symbolic links created with ln -s):\n";
	print $fh "    data = /data/DelGet/genomes\n";
	print $fh "#      The genome and gap files need to exist for each species, with the species_id from the newick tree .fa\n";
	print $fh "#      Genomes need to be spID.fa and gap files spID.gaps.tab\n";
	print $fh "#      Typically: hg38.fa and hg38.gaps.tab\n";
	print $fh "#      The gap files list the assembly gaps coordinates, in UCSC format (see the help)\n";
	print $fh "#      If no gap file available on UCSC for your assembly, just let the pipeline make them (or use the provided Util script fasta_get_gaps.pl)\n";
	print $fh "#      The genome index, total length files and gap files if needed will be written at the genome location: so it needs writing permissions\n";
	print $fh "\n";
	print $fh "# path of the folder that will contain all output files:\n";
	print $fh "    path = ./Del #This directory will be created if it does not exist\n";
	print $fh "\n";
	print $fh "# Set behavior for aligning:\n";
	print $fh "    still_align = yes #leave blank to skip that step\n"; 
	print $fh "\n";
	print $fh "# Set behavior for analyzing gaps:\n";
	print $fh "    analyze_gaps = yes #leave blank to skip that stepd\n"; 
	print $fh "\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# STEP OF GETTING REGIONS FROM OUTGROUP TO ANALYZE FOR DELETIONS: can be random or loaded from a file\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# Set an input file to load from file. Leave blank to get regions randomly.\n";
	print $fh "    file = \n"; 
	print $fh "#      Can only be a Repeat Masker output file (.out at the end mandatory) for now\n";
	print $fh "#      And has to be a masking of the outgroup, unless only 2 species are considered,\n";
	print $fh "#      in which case the 'outgroup' is the masked species. Rounds will be done randnb by randnb,\n";
	print $fh "#       or 1500 by 1500 if randnb is not set.\n";
#	print $fh "#      The format will be specified by the extension:\n";
#	print $fh "#         end by .out for a Repeat Masker output\n"; 
#	print $fh "#         end by .bed for a bed file, with at least 3 columns: chr \\t start \\t end \\t ID \\t...\n";
#	print $fh "#            (if a unique ID is defined it will be used, otherwise it will be made from coordinates)\n";
#	print $fh "#         end by .gff for a gff or gff3 file, where columns are: chr \\t XX \\t XX \\t start \\t end \\t ...\n"; 
#	print $fh "#         end by .gtf for a gtf file (will be imported like a gff file since the first columns are similar)\n"; 
	print $fh "# Set how anchors are defined (not relevant if regions are loaded randomly)\n";
	print $fh "    window = 1\n"; 
	print $fh "#      Set to 1 if anchors should flank each feature, such as:\n";
	print $fh "#         [TE1.anchor5'][==TE1==][TE1.anchor3'][=====//=====][TE2.anchor5'][==TE2==][TE2.anchor3']\n";
	print $fh "#      Set to 2 if they should flank consecutive features, such as:\n";
	print $fh "#         [anchor5'][==TE1==][=====//=====][==TE2==][anchor3']\n";
	print $fh "#The idea behind this is to see if non allelic recombination events can be detected,\n";
	print $fh "#         so the .out file needs to only contain the TEs of interest\n";
	print $fh "\n";
	print $fh "# length of anchors, in nt\n";
	print $fh "    anch_len = 100 #100 is a good number; 80 gave similar results\n";
	print $fh "\n";
	print $fh "# Highest ratio (second highest score / highest score). This will allow to filter for hits that are not high enough\n";
	print $fh "    maxratio = 0.9 #this means: second highest score / highest score < 0.9 (ie second score is max 90% of highest)\n";
	print $fh "  \n";
	print $fh "# Minimum length of a N strech for it to be considered as an assembly gap\n";
	print $fh "    mingaplen = 50 #assembly gaps are typically 100nt\n";
	print $fh "\n";
	print $fh "# [FYI - CHANGING HERE HAS NO IMPACT] this is how the maximum length of the anchor is set (+2nt)\n";
	print $fh "#   a_max_len = anch_len + 50*(anch_len)/100 # This is to avoid huge span of the hit of anchors in other genome\n";
	print $fh "  \n";
	print $fh "# [FYI - CHANGING HERE HAS NO IMPACT] this is how the total length of regions (+2nt) is calculated\n";
	print $fh "#   minlen = anch_dist + 2*anch_len\n";
	print $fh " \n";
	print $fh "----------------------------------------------------------\n";
	print $fh "# VARIABLES FOR RANDOMIZATION (skip besides randnb if file is set)\n";
	print $fh "----------------------------------------------------------\n";
	print $fh "# Total length of alignment (in nt) required to end the pipeline\n";
	print $fh "    randtot = 100000\n"; 
	print $fh "#      Note that rounds will be done randnb by randnb and printed as the pipeline progresses,\n";
	print $fh "#      which allows to obtain preliminary results. To do so, run the Util script: DelGet_cat-outputs.pl\n"; 
	print $fh "#      with the directory Del as argument, even if the runs are not done.\n";
	print $fh "#      randtot will be use to approximate the number of randomizations per chromosomes/scaffolds (randtot / anch_dist)\n";
	print $fh "\n";
	print $fh "# Number per run that need to be successful before going to extracting and aligning sequences\n";
	print $fh "    randnb = 50\n";
	print $fh "#      After a while, it is possible that very little or no anchors manage to get through the filters\n";
	print $fh "#      Running in loop allows to still get outputs and results, even if script has to be killed for that reason or any other\n";
	print $fh "#      Just use the script DelGet_cat-outputs.pl in utilities to gather the gap lengths from the outputs.\n";
	print $fh "#      This is also needed when a file is loaded, to avoid huge fasta file to blat (in this case the default is 1500)\n";
	print $fh "\n";
	print $fh "# Number of random positions that will be treated per round\n";
	print $fh "# (anchor sequences extracted from outgroup, blat against target genomes, checking steps etc).\n";
	print $fh "    a_per_round = 500\n";
	print $fh "#      This step is repeated inside script --1-- until randnb is reached. \n";
	print $fh "#      Blat step is a limitation in the script and requires to put genome in memory every time\n";
	print $fh "#      therefore, if species are closely related or if assembly is good, you can lower this number,\n";
	print $fh "#      but given the usual success rate even in primates, at least randnb x2 is advised.\n";
	print $fh "#      Do x5 to x10 if more than 3 species, if less good assembly or more distant species,\n"; 
	print $fh "#      or if they have a lot of recent TEs for example.\n";
	print $fh "\n";
	print $fh "# Length beetween anchors, in nt\n";
	print $fh "    anch_dist = 10000 #=distance in genome1 between the 2 anchors [anchor1]<----XXnt---->[anchor2]\n";
	print $fh "#      This will be use to approximate the number of randomizations per chromosomes/scaffolds (randtot / anch_dist)\n";
	print $fh "\n";
	print $fh "# Behavior for loading outputs generated during the run (to avoid overlaps):\n";
	print $fh "    if_OKregions = yes #Recommended\n";
	print $fh "#      Set to 'yes' or 'no':\n";
	print $fh "#      'yes' here combined with 'auto' below means that ALL previous OK regions will be concatenated by the pipeline.\n";
	print $fh "#      'yes' here and a path to a file below will load the file as the first previous output\n";
	print $fh "#      'no' here means that previous regions won't be checked for obverlaps (so the regions analyzed won't be independant!!)\n";
	print $fh "\n";
	print $fh "    OKregions = auto #Recommended\n";
	print $fh "#      Set to 'auto' or to a file:\n";
	print $fh "#      'auto' here means that previous regions that went through the filters will be pre-loaded to avoid overlaps,\n";
	print $fh "#         by concatenating Deletions.Run(-)1/_OKregions.tab and the Deletions.Run(-)2/_OKregions.tab (if any)\n";
	print $fh "#      If that is not desired, just move or rename the previous output directories (Deletion.RunX)\n";
	print $fh "#      Set to a specific file to append new regions to existing previous output (if there were some with different file names for example)\n";
	print $fh "\n";
	print $fh "----------------------------------------------------------\n";
	print $fh "# VARIABLES FOR EXTRACTION AND ALIGNMENT\n";
	print $fh "----------------------------------------------------------\n";
	print $fh "# Muscle does not like very long regions - will switch to Kalign if one sequence is > max_muscle\n";
	print $fh "    max_muscle = 30000 #Recommended\n";
	print $fh "\n";	
	print $fh "# [FYI - CHANGING HERE HAS NO IMPACT] this how the max length of sequence to extract is determined,\n";
	print $fh "#   to avoid aligning regions that are too long. Note: will try muscle if < max_muscle and then kalign if it crashes.\n";
	print $fh "#   multip = 3                                             #ex. anch_dist = 25kb => 75kb max\n";	
	print $fh "#   multip = 4 if anch_dist <= 10000                       #ex. anch_dist = 10kb => 40kb max\n";
	print $fh "\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# PATHS TO SOFTWARES\n";
	print $fh "#------------------------------------------------------------------------------------------------------------------\n";
	print $fh "# Blat stand alone, see http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3\n";
	print $fh "    blat = /home/software/ucsc/blat/blat\n";
	print $fh "\n";
	print $fh "# Muscle software, see http://www.drive5.com/muscle/downloads.htm \n";
	print $fh "    muscle = /home/software/muscle3.8.31/muscle3.8.31\n";
	print $fh "#      When this one is use: all sequences are re-written in upercase by default.\n";
	print $fh "\n";	
	print $fh "# Kalign software [needed only for large alignments], see http://www.biomedcentral.com/1471-2105/6/298/\n";
	print $fh "    kalign = /home/software/Kalign2/kalign\n";
	print $fh "#      This software keeps the original case.\n";
	print $fh "\n";
	print $fh "# [FYI - CHANGING HERE HAS NO IMPACT]: alignment software will be decided based on anch_dist\n";
	print $fh "#   alnsoft = muscle if anch_dist < 30000 nt\n";
	print $fh "#   alnsoft = kalign if anch_dist >= 30000 nt\n";
	print $fh "\n";
	print $fh "# BEDtools bin directory, see https://github.com/arq5x/bedtools2\n";
	print $fh "    bedtools = /home/software/bedtools2/bin\n";
	print $fh "#      Notes on BEDtools command lines:\n";
	print $fh "#         - Two gaps between two species are considered to be specific if they overlap less than 80% for gaps 1-30nt and 85% for >30nt\n";
	print $fh "#           e.g. option -f 0.80 of intersecBed = flexibility of 1nt per 10nt of gaps\n";
	print $fh "#              These 2 gaps would be shared (17nt of overlap is > 80% of the gap)\n";
	print $fh "#                 specie 1 ATGCATGCATGCATGCA--------------------ATGCATGCATGCATGCATGCATGC\n";
	print $fh "#                 specie 2 ATGCATGCATGCAT--------------------GCAATGCATGCATGCATGCATGCATGC\n";
	print $fh "#              But these 2 would be specific (15nt of overlap is < 80% of the gap)\n";
	print $fh "#                 specie 1 ATGCATGCATGCATGCA--------------------ATGCATGCATGCATGCATGCATGC\n";
	print $fh "#                 specie 2 ATGCATGCATGCAT------------------ATGCAATGCATGCATGCATGCATGCATGC\n";
	print $fh "#         - This minimum overlap is reciprocal\n";
	print $fh "#           (if BIG gap in one specie, and a small one in the other, likely not from the same deletion event)\n";
	print $fh "#           => option -r of intersectBed\n";
	print $fh "#              To be able to use -r, intersectBed -v instead of subtractBed is used\n";
	print $fh "#              (-v means keep only stuff that don't overlap)\n";
	print $fh "\n";          
# 	print $fh "# Path to the Repeat Masker software, installed and configured (leave blank for no masking)\n";
# 	print $fh "    RMSOFT = /home/software/RepeatMasker/RepeatMasker\n";
# 	print $fh "\n";
# 	print $fh "----------------------------------------------------------\n";
# 	print $fh "# REPEAT MASKING\n";
# 	print $fh "----------------------------------------------------------\n";
# 	print $fh "# Related to Repeat Masking:\n";
# 	print $fh "# Librarie(s):\n";
# 	print $fh "    mask = #list of fasta files (leave blank for no masking)\n";
# 	print $fh "#          to mask with several libraries, separate them by commas and in between parentheses, in the mask and the masking_folders list\n";
# 	print $fh "#          Note that -lib is needed in running repeat masker => use Repeat Masker script in their utilities\n";
# 	print $fh "#          to convert their embl to fasta if you would like to mask with the whole library\n";
# 	print $fh "# Output directories for the masked sequences\n";
# 	print $fh "    masking_folders = \n";
# 	print $fh "#          as many masking_folders elements as library files are required\n";
# 	print $fh "\n";
	close $fh;	
	print STDERR "\n   Configuration file example printed in $CONF - once edited, run:\n";
	print STDERR "        perl DelGet.pl $CONF\n\n";	
	return 1;
}

#------------------------------------------------------------------------------
sub load_config {
	open(my $fh, "<", $CONFIG_FILE) or confess "\t    ERROR: SUB print_config: could not open to write $CONFIG_FILE!\n";
	while(defined(my $line = <$fh>)) {	
		chomp $line;
		$line =~ s/\s//g;
		next if ((substr($line,0,1) eq "#") || ($line !~ /\w/) || ($line !~ /=/));
		my @line = split("=",$line);
		$line[1] = "" unless ($line[1]);
		$line[1] =~ s/#.*$//;
		$CF{$line[0]}=$line[1];
	}
	close($fh);
	
	#Now the ones that are calculated
	$CF{'a_max_len'} = $CF{'anch_len'} + (50 * $CF{'anch_len'} /100);
	if ($CF{'file'}) {
		$CF{'alnsoft'} = $CF{'kalign'}; #safer that way
		$CF{'r_max_len'} = 100000; #100 kb is... long enough.
		$CF{'randtot'} = `wc -l < $CF{'file'}`;		
		chomp($CF{'randtot'});
		$CF{'randnb'} = 1500 if (! $CF{'randnb'});
	} else {
		if ($CF{'anch_dist'} < $CF{'max_muscle'}) {
			$CF{'alnsoft'} = $CF{'muscle'};
		} else {
			$CF{'alnsoft'} = $CF{'kalign'};    
		}
		$CF{'randtot_c'} = $CF{'randtot'} / $CF{'anch_dist'};
		$CF{'minlen'} = $CF{'anch_dist'} + (2 * $CF{'anch_len'});
		$CF{'multip'} = 3;
		$CF{'multip'} = 4 if ($CF{'anch_dist'} <= 10000); #10kb intially => 40kb max\n";
		$CF{'r_max_len'} = $CF{'anch_dist'} * $CF{'multip'};
	}	
	return 1;
}

#------------------------------------------------------------------------------
sub load_newick_tree {
	if (-f $CF{'tree'}) {
		print $MFH " LOADING TREE from $CF{'tree'}\n";
		my $nwck = `cat $CF{'tree'}`;
		chomp($nwck);
		$nwck =~ s/\n//;
		$CF{'tree'} = $nwck;
	} 	
	print $MFH " --- Species tree is as follow:\n";
	print $MFH "     $CF{'tree'}\n";
	return 1;
}

#----------------------------------------------------------------------------
sub parse_newick_tree {	
	if ($CF{'tree'} =~ /\(([\w\-_,]+?)\)/) {
		$NWCK{$LEAF} = $CF{'tree'};
		$NWCK{$LEAF} =~ s/^.*\(([\w\-_,]+?)\).*$/$1/;				
		get_spids(); #load list of species IDs of this leaf		
		$CF{'tree'} =~ s/(^.*)\(([\w\-_,]+?)\)(.*$)/$1$LEAF$3/;		
		$MAX = $LEAF;
		$LEAF++;		
	} else {
		$ITER = 1;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_spids {
	if ($NWCK{$LEAF} =~ /,/) {
		my @ids = split(",",$NWCK{$LEAF});
		foreach my $id (@ids) {
			if ($SPIDS{$id}) {
				my @id = @{$SPIDS{$id}};
				push (@{$SPIDS{$LEAF}},@id);
			} else {
				push (@{$SPIDS{$LEAF}},$id);
			}
		}
	} else {
		my $id = $NWCK{$LEAF};
		$id = $SPIDS{$id} if ($SPIDS{$id});
		$SPIDS{$LEAF}=$id;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_groups_from_ids {
	while (my($key, $value) = each %SPIDS) {
		$GROUP{@{$value}}=$key;
		my @val = @{$value};		
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_group_parents {	
	$PAR{$MAX}="root"; #max lvl = all species
	foreach my $lvl (sort { $a <=> $b } keys %NWCK) { 
		my @sp = @{$SPIDS{$lvl}};
		#parent is the last group that contained the species ids in its list
		for (my $i = $lvl+1; $i < $MAX+1; $i++) {
			my @thsp = @{$SPIDS{$i}};
			my $found;
			$found = check_array(\@sp,\@thsp);
			$PAR{$lvl}=$GROUP{@thsp} if ($found);
			last if ($found);
		}
	}
	return 1;
}

#----------------------------------------------------------------------------
sub check_array {
    my ( $test, $source ) = @_;
    my %exists = map { $_ => 1 } @$source;
    foreach my $ts ( @{$test} ) {
        return if ! $exists{$ts};
    }
    return 1;
}

#----------------------------------------------------------------------------
sub get_genome_files_list {
	foreach my $sp (@SPIDS) {
		my $fa = "$CF{'data'}/$sp.fa";
		if (-e $fa) {
			$DATA{$sp}{'f'}=$fa;
			$DATA{$sp}{'n'}="$sp.fa";
		} else {		
			die "\t    ERROR - genome file for $sp could not be found - $fa does not exist?\n\n";
		}
# 		if ($CF{'still_align'}) {
# 			my $fa2 = "$CF{'data'}/$sp.2.fa";
# 			`ln -s $CF{'data'}/$sp.fa $fa2` if (! -e $fa2);
# 			$DATA{$sp}{'e'}=$fa2;
# 		}
	}
	return 1;
}

#------------------------------------------------------------------------------
sub get_fa_infos {
	my $fa = $DATA{$OG}{'f'};
	my $minlen = $CF{'minlen'}*2; #skip sequences that are too short => set to x2 the region length
	my $lout = $fa.".tot-$CF{'randtot_c'}.min-$minlen.infos.txt";  
	if (! -e $lout) {
		print $MFH "        Getting sequence lengths > $minlen ($CF{'minlen'}*2) and printing data in $lout\n";
		my %tmp = ();
		my $id = "";
		my $ln = 0;		
		my $c = 0;			
		open (my $fh, "<", $fa) or confess "\t    ERROR - could not open to read file $fa $!\n\n";
		while (defined(my $l = <$fh>)) {
			chomp ($l);
			if (substr($l,0,1) eq ">") {
				#check length, if not first seq
				if ($c == 1 && $ln > $minlen) { 
					$DATA{$OG}{'l'}+=$ln;
					$tmp{$id}=$ln;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$l);
				$id = $id[0];
				$id =~ s/>//;
				$ln = 0;
			} else {
				$ln+=length($l);
			}
		}
		#deal with last seq
		if ($ln > $minlen) {
			$DATA{$OG}{'l'}+=$ln;
			$tmp{$id}=$ln;
		}
		close $fh;
		
		print $MFH "NOW PRINTING\n";
		
		#Now print the data and get the infos
		open (my $lfh, ">", $lout) or confess "\t    ERROR - could not open to write file $lout $!\n\n";
		foreach my $id (keys %tmp) {
			my $ratio = int($CF{'randtot_c'} / $DATA{$OG}{'l'} * $tmp{$id});			
			$ratio = 1 if ($ratio == 0); #bug fix 2015 05 13
			print $lfh "$id\t$tmp{$id}\t$ratio\n";
			my @infos = ($tmp{$id},$ratio);
			$INFOS{$id} = \@infos;
			push(@DBIDS,$id);
		}
		close $lfh;
		%tmp = ();
	} else {
		print $MFH "        Loading data from $lout\n";
		open (my $lfh, "<", $lout) or confess "\t    ERROR - could not open to read $lout $!\n";
		while (defined(my $l = <$lfh>)) {
			chomp ($l);
			my ($id,$len,$ratio) = split (/\t/,$l);
			my @infos = ($len,$ratio);
			$INFOS{$id} = \@infos;
			push(@DBIDS,$id);
		}
		close $lfh;
	}
	return 1;
}

#----------------------------------------------------------------------------
sub get_gap_files_list {
	foreach my $sp (@SPIDS) {
		my $fa = "$CF{'data'}/$sp.fa";
		my $gap = `find $CF{'data'}/$sp*.gaps*tab`;
		chomp ($gap);
		if (-e $gap) {	
			$DATA{$sp}{'g'} = $gap;
		} else {
			print $MFH "   Gap file for $sp could not be found - $CF{'data'}/$sp*.gaps*tab does not exist\n";
			print $MFH "     => obtain, with command line:\n";
			print $MFH "        perl $BIN/Utils/fasta_get_gaps.pl -i $CF{'data'}/$fa -l 50\n";
			`perl $BIN/Utils/fasta_get_gaps.pl -i $CF{'data'}/$fa -l 50`;
			$DATA{$sp}{'g'}="$fa.gaps.min50.tab";
		}
	}
	return 1;
}

#------------------------------------------------------------------------------
sub load_gaps {
	foreach my $sp (@SPIDS) {
		open(my $fh, "<", $DATA{$sp}{'g'}) or confess "\t    ERROR: could not open to read $DATA{$sp}{'g'}!\n";
		#0		1							2			3			4	5	6
		#bin	chrom						chromStart	chromEnd	ix	n	size	type	bridge
		#.		gi|432120764|gb|KB097841.1|	268			451			.	.	184		.		.
		while(defined(my $l = <$fh>)) {	
			chomp $l;
			next if ($l =~ /#bin/ || $l !~ /\w/);
			my @gl = split(/\t/,$l);
			$GAPS{$sp}{$gl[1]} .= "$gl[2],$gl[3],";
		}
		close $fh;
	}	
	return 1; 
}

#------------------------------------------------------------------------------
sub load_file {
	my $file = $CF{'file'};
	my $len = $CF{'anch_len'};
	my $w = $CF{'window'};
	my @list = ();
	open(my $fh, "<", $file) or confess "\t    ERROR: SUB load_coords: could not open to write $file!\n";
	while(defined(my $l = <$fh>)) {	
		chomp $l;
		my @l = ();
		$l =~ s/^\s+//;
		next if ($l =~ /^[Ss]core|^SW|^#/ || $l !~ /\w/);
		@l = split(/\s+/,$l);
		my @coords = ($l[4],$l[5],$l[6],$l[9]);			
		push(@list,\@coords);	
	}
	close $fh;
	
	#now deal with the whole window thing - I only need the anchors' coordinates
	$w--; #so that same window is a 0
	my $j = 1;
	my $c = 1;
	my %te_per_chr = ();
	for (my $i=0;$i<=$#list-$w;$i++) {
		my $chr = $list[$i]->[0];
		if (exists $te_per_chr{$chr}) {
			$te_per_chr{$chr}++;
		} else {
			$te_per_chr{$chr}=0;
		}				
		my $en5 = $list[$i]->[1]; #end of 5' anchor
		my $st5 = $en5 - $len+1; #start of the 5' anchor, len is really $len 
		my $st3 = $list[$i+$w]->[2];  #start of 3' anchor; it will be the same feature if w = 1
		my $en3 = $st3 + $len-1; #end of 3' anchor; it will be the same feature if w = 1
		$COORDS{$c}{$chr}{$te_per_chr{$chr}}[0] = $st5;
		$COORDS{$c}{$chr}{$te_per_chr{$chr}}[1] = $en5;
		$COORDS{$c}{$chr}{$te_per_chr{$chr}}[2] = $st3;
		$COORDS{$c}{$chr}{$te_per_chr{$chr}}[3] = $en3;
		$COORDS{$c}{$chr}{$te_per_chr{$chr}}[4] = "$list[$i]->[3].$list[$i+$w]->[3]";	
		if ($j == $CF{'randnb'}) {
			$c++;
		}
		$j++;
	}
	return 1;
}

#------------------------------------------------------------------------------
sub make_out_dirs {
	until (!-d "$CF{'path'}/Deletions.$C"){
		$C++;
	}
	$CURRDEL = "$CF{'path'}/Deletions.$C";
	mkdir ($CURRDEL, 0755) or confess "\t    ERROR: can not mkdir $CURRDEL $!";	
	
	$ANCHORS = "$CURRDEL/data_anchors-posi-fa-blatout";
	mkdir ($ANCHORS, 0755) or confess "\t    ERROR - Couldn't mkdir $ANCHORS $!";

	$BLATS = "$CURRDEL/data_highest-scores";		
	mkdir ($BLATS, 0755) or confess "\t    ERROR - Couldn't mkdir $BLATS $!";
	
	$ALIGN = "$CURRDEL/_ExtractAlign";
	mkdir ($ALIGN, 0755) or confess "\t    ERROR - Couldn't mkdir $ALIGN $!";
	
	$GAPS = "$CURRDEL/_GapAnalysis";
	mkdir ($GAPS, 0755) or confess "\t    ERROR - Couldn't mkdir $GAPS $!";
	
	return 1;
}



#------------------------------------------------------------------------------
# GET & PROCESS REGIONS
#------------------------------------------------------------------------------
sub get_previous_OKreg {
	# Structure of _OKregions.tab output file starts like this:
	# 0		1		2		3		4		5		6		7		8		9
	# ID	type	Gname	Gstart	Gend	type	Gname	Gstart	Gend	DIST
	my $name = "_OKregions.tab";
	my $nameall = "_OKregions.cat.tab";
	my $toload;
	#Nothing to do with $C == 0
	#So start with 1
	if ($C != 0) {
		if ($C == 1) { 
			if ($CF{'OKregions'} eq "auto") {
				#load regions from Deletions.0
				$toload = "$CF{'path'}/Deletions.0/$name";
			} else { 
				#need to cat the defined regions in config file with OKregions in Deletions.0
				$toload = "$CF{'path'}/Deletions.0/$nameall";
				system "cat $CF{'path'}/Deletions.0/$name $CF{'OKregions'}> $CF{'path'}/Deletions.0/$nameall";			
			}
		} else { #more than 1 => cat the -2 nameall and new one
			my $c_prev1 = $C - 1;
			my $c_prev2 = $C - 2;
			$toload = "$CF{'path'}/Deletions.$c_prev1/$nameall";
			my $prev2;
			if ($c_prev2 == 0) {
				$prev2 = "$CF{'path'}/Deletions.$c_prev2/$name"; #no cat file in Deletions.0
			} else {
				$prev2 = "$CF{'path'}/Deletions.$c_prev2/$nameall";
			}
			system "cat $prev2 $CF{'path'}/Deletions.$c_prev1/$name > $toload";
		}
		print $MFH "       => $toload will be loaded to avoid overlaps\n";
	
		#Now, load OK region file unless not relevant, ie no previous runs
		open(my $pfh, ">", $PREVREG) or confess "\t    ERROR - can not open to write $PREVREG $!";
		print $pfh "#chr\tend\n\n";
		open(my $okfh, "<", $toload) or confess "\t    ERROR - can not open to read previous _OKregions output file $toload $!";
		my $reg_r;
		while(defined(my $l = <$okfh>)) {	
			chomp $l;
			next if ($l =~ /^#/ || $l !~ /\w/);
			my @l = split(/\t/,$l);
			print $pfh "$l[2]\t$l[10]\n";			
			# get the round to start with
			my $curr_r = $l[0];
			$curr_r =~ s/^reg(.*)-.*$/$1/; #extract this line round number
			$reg_r = $curr_r unless ($reg_r); #first round, has to be remembered
			#memorize region ID, if larger (that way still OK even if not numerical order in the file, i.e. excel sorts by alphabetical)
			$reg_r = $curr_r if ($reg_r < $curr_r); 
		}
		$R = $reg_r + 1; #reinitialize $r to this last round number +1 => still unique IDs for regions
		print $MFH "          -> last iteration was $reg_r => new regions will start at $R\n";
		close $pfh;
		close $okfh;
	}
	return 1;
}

#------------------------------------------------------------------------------
sub prep_curr_out_files {
	open(my $ffh,">",$OKREG) or die "\t    ERROR - can not open to write $OKREG $!";
	print $ffh "#REGION\t";
	foreach my $sp (@SPIDS) {
	    print $ffh "$sp\t\t\t\t\t\t\t\t\t\t\t\t\t";
	}
	print $ffh "\n";
	print $ffh "#region_ID\t";
	foreach my $sp (@SPIDS) {
	    print $ffh "type\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\tstrand\tDIST\t";
	    #Will just have three "." values for the outgroup, it's fine
	}
	print $ffh "\n";
	close $ffh;
	
	if ($CF{'still_align'}) {
		#POSI output:
		#region1     scaffold/chr     start     end     +     Mluc     /mnt/disk3/genomes/Myotis_lucifugus/Myotis_7x.fasta
		open(my $pfh,">",$POSI) or confess "\t    ERROR - can not open to write $POSI $!";
		print $pfh "#region\tchr\tstart\tend\tstrand\tgenome_ID\tgenome_location\n\n";
		close $pfh;
	}	
	return 1;
}

#------------------------------------------------------------------------------
sub load_previous_OKreg {
	print $RFH " --- Loading previous OK regions $PREVREG\n";
	open(my $fh, "<", $PREVREG) or confess "\t    ERROR - can not open to read file $PREVREG $!";
	while(defined(my $l = <$fh>)) {	
		chomp $l;
		next if ($l =~ /^#/ || $l !~ /\w/);
		my @l = split(/\t/,$l);
		if (exists $ALREADYRAND{$l[0]}) {
			$ALREADYRAND{$l[0]} .= ",$l[1]";
		} else {
			$ALREADYRAND{$l[0]} = $l[1];
		}
	}
	close $fh;
	return 1;
}	

#------------------------------------------------------------------------------
sub get_regions_from_file {
	print $RFH " --- Looping through and write positions\n";
	my $i = 0;
	foreach my $chr (keys %{$COORDS{$R}}) {
		my $coordsnb = scalar(keys %{$COORDS{$R}{$chr}});		
		for (my $i=0;$i<$coordsnb;$i++) {		
			#Check for assembly gaps; should not happen but could, if too close to one.
			#Sub takes 3'end first and then 5'end
			my $ifnext = check_overlap_gap($GAPS{$OG}{$chr},$COORDS{$R}{$chr}{$i}[3],$COORDS{$R}{$chr}{$i}[0]);
			if ($ifnext eq "yes") {
				$FAILED_ANCHOR++ ;
				next;
			}	
			
			#OK => print in posifiles 
			open(my $pfh,">>",$DATA{$OG}{'posi'}) or confess "\t    ERROR - can not open to write $DATA{$OG}{'posi'} $!";
			my $regname = "reg.$chr.$COORDS{$R}{$chr}{$i}[1].$COORDS{$R}{$chr}{$i}[2]._.$COORDS{$R}{$chr}{$i}[4]";			
			print $pfh $regname."#5\t$chr\t$COORDS{$R}{$chr}{$i}[0]\t$COORDS{$R}{$chr}{$i}[1]\n";
			print $pfh $regname."#3\t$chr\t$COORDS{$R}{$chr}{$i}[2]\t$COORDS{$R}{$chr}{$i}[3]\n";
			$i++;
			close $pfh;
		}
	}	
	print $RFH "     => $FAILED_ANCHOR regions skipped because overlapping asssembly gaps\n";
	return 1;
}

#-------------------------------------------------------------------------------
sub get_regions {
	print $RFH " --- Getting $CF{'a_per_round'} random positions\n";
# anchors from the outgroup:
# 	$DATA{$OG}{'posi'} = "$ANCHORS/$R.posi.$OG.tab";
# 	$DATA{$OG}{'anchors'} = "$ANCHORS/$R.anchors.$OG.fa";
	my %infos = %INFOS;
	my $i = 0;
	my %rand_ends = (); #save what is already randomized in here to avoid overlaps (for next rounds it is checked with printed regions)
	until ($i == $CF{'a_per_round'}) {
		#Randomize using list of IDs that went through length filter
		my $rdmIDposi = int(rand($#DBIDS)); 
		my $rdmID = $DBIDS[$rdmIDposi];
		my $ratio = $infos{$rdmID}->[1]; #corresponds to nb of random positions to do here, weightened on length
		my $size = $infos{$rdmID}->[0] - ($CF{'minlen'}-2); #-2 bc +1 before, and bc minlen will be added to all values	
		if ($ratio < 1 || $size < 1) { #shouldn't happen since now IDs are only the ones of contigs larger than $minlen, but just in case keep this
			$FAILED_ANCHOR++;
			#print $RFH "ratio_or_size<1=>next\n";
			next;
		}	
		my $three_end = int(rand($size)) + 1 + ($CF{'minlen'}-2); #get in base1

		# get other coordinates and store the end
		my $three_start = $three_end - $CF{'anch_len'} + 1; #len of anchor = 100nt
		my $five_end = $three_start - $CF{'anch_dist'} + 1; #len in between is really 300nt
		my $five_start = $five_end - $CF{'anch_len'} + 1; #len of anchor = 100nt
		if ($rand_ends{$rdmID}) {
			$rand_ends{$rdmID} .= ",$three_end";
		} else {
			$rand_ends{$rdmID} = $three_end;
		}	
				
		#3) avoid overlap with any already picked region
		my $ifnextrand = "no";
		# - first, do not overlap current randomizations
		$ifnextrand = check_overlap_reg($rand_ends{$rdmID},$three_end,$five_start) if (exists $rand_ends{$rdmID});
		if ($ifnextrand eq "yes") {
			$FAILED_ANCHOR++;
			next;
		}				
		# - second, check regions previously stored
		($ifnextrand) = check_overlap_reg($ALREADYRAND{$rdmID},$three_end,$five_start) if (exists $ALREADYRAND{$rdmID});
		if ($ifnextrand eq "yes") {
			$FAILED_ANCHOR++;
			next;
		}
		# => not overlapping with previous stuff => OK to continue 
		
		#4) Checking assembly gaps; start and end storred this time b/c size of gaps may change		
		($ifnextrand) = check_overlap_gap($GAPS{$OG}{$rdmID},$three_end,$five_start) if (exists $GAPS{$OG}{$rdmID});
		if ($ifnextrand eq "yes") {
			$FAILED_ANCHOR++;
			next;
		}
				
		#everything OK => print in posifile 
		print $RFH "       => $i regions, OK\n" if ($i >= $CF{'a_per_round'});		
		open(my $pfh,">>",$DATA{$OG}{'posi'}) or confess "\t    ERROR - can not open to write $DATA{$OG}{'posi'} $!";
		print $pfh "reg".$R."-".$i."#5"."\t$rdmID\t$five_start\t$five_end\n";
		print $pfh "reg".$R."-".$i."#3"."\t$rdmID\t$three_start\t$three_end\n";
		$i++;
		close $pfh;
		#since there was 1 random succesful for this one, remove 1 #Moved here bug fix 2016 07 27
		$infos{$rdmID}->[1]--;
	}
	print $RFH "     => failed randomizations = $FAILED_ANCHOR\n";
	return 1;	
}	

#------------------------------------------------------------------------------
sub check_overlap_reg {
	my ($already,$three_end,$five_start) = @_;
	my $ifnext = "no";
	my @already = split(",",$already);
	foreach my $prev_end (@already) {
		my $prev_start = $prev_end - (2*$CF{'anch_len'}) - $CF{'anch_dist'} + 3;
		if (($three_end < $prev_end && $three_end > $prev_start) || ($five_start < $prev_end && $five_start > $prev_start)) {
			$ifnext = "yes"; #there is overlap, no need to check the rest => will skip this region
			last;
		}
	}
	return ($ifnext);
}

#------------------------------------------------------------------------------
sub check_overlap_gap {
	my ($storedgaps,$three_end,$five_start) = @_;
	my $ifnext = "no";
	return ($ifnext) unless ($storedgaps);
	my @gaps = split(",",$storedgaps);
	for (my $g=0;$g<=$#gaps;$g=$g+2) {
		my $gap_start = $gaps[$g];
		my $gap_end = $gaps[$g+1];
		if ($gap_end > $five_start && $gap_start < $three_end) {
			$ifnext = "yes"; #there is overlap, no need to check the rest => will skip this region
			last;
		}	
	}
	return ($ifnext);
}

#------------------------------------------------------------------------------
sub extract_anchor_seq {	
#outgroup files:
# 	$DATA{$OG}{'f'}
# 	$DATA{$OG}{'posi'}
# 	$DATA{$OG}{'anchors'}
	print $RFH " --- Extracting fasta sequences of the anchors (from $OG)\n";
	my $db = Bio::DB::Fasta->new($DATA{$OG}{'f'}) 
	         or confess "\t    ERROR - could not create Bio::DB::Fasta object from $DATA{$OG}{'f'} $!\n";
	print $RFH " --- Extracting sequences...\n";
	
	#extract
	my $og_infos = ();
	open(my $fh,"<", $DATA{$OG}{'posi'}) or confess "\t    ERROR - can not open file $DATA{$OG}{'posi'} $!";
	while(defined(my $l = <$fh>)) {	
		next unless ($l =~ /\w/);	#avoid blank lines
		chomp($l);
		my ($rname,$Gname,$Gstart,$Gend)= split(/\t/,$l);	
		my $subSeq = $db->seq($Gname,$Gstart,$Gend) or print $RFH "        WARN: $Gname was not found in $DATA{$OG}{'f'}\n";	
		my $seqobj = Bio::Seq->new( -display_id => $rname, 
									-seq        => $subSeq);
		my $aseqio = Bio::SeqIO->newFh(-format => 'Fasta', 
								       -file=>">>$DATA{$OG}{'anchors'}") 
				     or confess "Failed to create SeqIO FH object from $DATA{$OG}{'anchors'} $!\n";
	
		print $aseqio $seqobj;
	
		#keep coords in memory for final output
		my ($name,$type) = split("#",$rname);	
		$OGINFOS{$name}{$type} = "$Gname\t$Gstart\t$Gend";	
	}
	close $fh;	
	return 1;
}

#------------------------------------------------------------------------------
sub blat_anchors {
#outgroup files:
# 	$DATA{$OG}{'f'}
# 	$DATA{$OG}{'anchors'}
#for other species:
# 	$DATA{$sp}{'blat'}
	print $RFH " --- Blat these anchors against other genomes\n";
	foreach my $sp (@SPIDS) {
		next if ($sp eq $OG);
		print $RFH "       -> $sp\n";
		$DATA{$sp}{'blat'} = "$DATA{$OG}{'anchors'}.blat.$sp";
#		print $RFH "          $CF{'blat'} $DATA{$sp}{'f'} $DATA{$OG}{'anchors'} -minScore=10 -minIdentity=10 $DATA{$sp}{'blat'}\n";
		print $RFH "          $CF{'blat'} $DATA{$sp}{'f'} $DATA{$OG}{'anchors'} -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 $DATA{$sp}{'blat'}\n";
#		`$CF{'blat'} $DATA{$sp}{'f'} $DATA{$OG}{'anchors'} -minScore=10 -minIdentity=10 $DATA{$sp}{'blat'} 2>> $RLOG`;
		`$CF{'blat'} $DATA{$sp}{'f'} $DATA{$OG}{'anchors'} -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 $DATA{$sp}{'blat'} 2>> $RLOG`;
	}	
	return 1;
#STRUCTURE OF PSL BLAT OUTPUT AND VARIABLE NAMES TO CALL VALUES
#FROM http://doc.bioperl.org/releases/bioperl-1.4/Bio/SearchIO/psl.html
#0			1			2				3		4				5				6				7				8			9		10			11			12		13		 14			15			16		17	
#$matches  $mismatches  $rep_matches  $n_count  $q_num_insert  $q_base_insert  $t_num_insert   $t_base_insert   $strand   $q_name   $q_length   $q_start  $q_end   $t_name   $t_length  $t_start   $t_end   $block_count  
#18				 19			 20
#$block_sizes    $q_starts   $t_starts

#0		1		2		3	4		5		6		7		8		9			10		11		12	13			14		15		16	17		18			19		 20
# match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	# match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count			
#--------------------------------------------------------------------------------------------------------------------------------------------------------------																	
# 17	0	0	0	0	0	0	0	+	GL429775_6626316_6626416_300bp_seq2_luci5'	101	83	100	gi|432119242|gb|KB098985.1|	2558592	1429419	1429436	1	17,	83,	1429419,
}

#------------------------------------------------------------------------------
sub parse_blat {
#File to parse:
# 	$DATA{$sp}{'blat'}
	print $RFH " --- Parsing the blat outputs\n";
	foreach my $sp (@SPIDS) {
		next if ($sp eq $OG);
		print $RFH "       -> $sp\n";

		#going through input + writing output
		my %highest = ();
		my %howhigh = ();
		my %blat = ();
		my $lc = 0;
		open(my $fh,"<",$DATA{$sp}{'blat'}) or confess "\t    ERROR - can not open blat output file $DATA{$sp}{'blat'} $!";
		while(defined(my $l = <$fh>)) {	
			chomp $l;
			unless ($lc < 5) {
				my @l = split(/\s+/,$l);
				my ($matches,$mismatches,$rep_matches,$q_name,$q_length) = ($l[0],$l[1],$l[2],$l[9],$l[10]);
				my $sc = sprintf "%.2f", ( 100 * ( $matches + $mismatches + $rep_matches ) / $q_length );
				my ($reg,$type) = split("#",$q_name);
				my $id = $reg."#".$sp."#".$type;
				if (exists $highest{$id}) {
					if ($sc > $highest{$id}) {
						$howhigh{$id} = $highest{$id} / $sc;
						$highest{$id} = $sc;
						$blat{$id} = $l;
					}
				} else { #first line basically
					$highest{$id} = $sc;
					$howhigh{$id} = 0;
					$blat{$id} = $l;
				}
			}
			$lc ++;
		}
		close $fh;
		#write an output to keep track, and keep file names
		my @HSfiles_temp = ();
		foreach my $scID (sort keys %highest) {
			my ($reg,$sp,$type) = split("#",$scID);
			my $HSonly = "$BLATS/$reg.HighestScores.tab";	
			open(my $HSonly_fh,">>",$HSonly) or confess "\t    ERROR - can not open to write $HSonly $!";		
			print $HSonly_fh "$blat{$scID}\t$reg\t$type\t$sp\tscore=\t$highest{$scID}\tratio_with_2nd=\t$howhigh{$scID}\n";	
			close $HSonly_fh;
			#Save file name
			push (@HSfiles_temp,$HSonly);	
		}
		push (@HSFILES,@HSfiles_temp);	
	}
	return 1;
	#STRUCTURE OF HighestScores.tab
	#0	1	2	3	4	5	6	7	8	9		10	11	12	13							14		15		16		17	18		19	20			21		22	23	24		25	26				27
	#97	4	0	0	0	0	0	0	+	reg2#3	101	0	101	gi|432108408|gb|KB104502.1|	6866958	2472779	2472880	1	101,	0,	2472779,	reg2	3	sp	score=	100	ratio_with_2nd=	0.2277
}

#------------------------------------------------------------------------------
sub check_and_print {
	my $ok = 0;
	print $RFH " --- Check each regions & print if passes the filters\n";
	print $RFH "     e.g. if hit, if same target, if score is high enough, if same strand, if overlap gaps in assembly...\n";
	open(my $ffh,">>",$OKREG) or confess "\t    ERROR - can not open to write $OKREG $!";
	open(my $pfh,">>",$POSI) or confess "\t    ERROR - can not open to write $POSI $!" if ($CF{'still_align'});
	for (my $i = 0; $i<=$#HSFILES; $i ++) {
		my $f = $HSFILES[$i];
		my %gen;
		my $reg;
		print $RFH "       -> $f\n";
		open (my $fh,"<",$f) or confess "\t    ERROR - can not open to read file $f $!";
		while(defined(my $l = <$fh>)) {	
			chomp($l);
			my @l = split(/\s+/,$l);
			my $sp = $l[23];
			my $type = $l[22];	                                       #  0		  1		   2	  3		  4	  5		 6
			my @vals = ($l[13],$l[15],$l[16],$l[8],$sp,$l[27],$l[25]); #($t_name,$t_start,$t_end,$strand,$sp,$ratio,$score)
			$reg = $l[21];
			$gen{$sp}{$type} = \@vals;
		}
		close $fh;
		
		#Check steps of the hits
		my ($ifnext,$failed) = check_anchor_hits(\%gen);
		if ($ifnext eq "yes") {
			print $RFH "            WARN: $reg failed because of: $failed\n";
			$FAILED_CHECKS++;
			next;
		}	
		$OK++;
		$ok++;
		
		#OK - print now
		print $RFH "            $reg went through filters => printing\n";	
		#print the region ID:	
		print $ffh "$reg\t";
		#then, for each sp: type\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\ttype\tGname\tGstart\tGend\tscore\tratio_with_next_highest_score\tstrand\tDIST				
		my ($regst,$regen,$dist);
		foreach my $sp (@SPIDS) { #outgroup is the first species of the list
			my $out = $CURRDEL."/_$sp.dist.tab";
			open(my $fh,">>",$out) or confess "\t    ERROR - can not open to write $out $!";	
			if ($sp eq $OG) {
				#$OGINFOS{$name}{$type} = "$Gname\t$Gstart\t$Gend";
				my @info5 = split(/\t/,$OGINFOS{$reg}{'5'});
				my @info3 = split(/\t/,$OGINFOS{$reg}{'3'});
				$regst = $info5[1];
				$regen = $info3[2];
				$dist = $regen - $regst +1;
				print $ffh "5\t$OGINFOS{$reg}{'5'}\tna\tna\t3\t\t$OGINFOS{$reg}{'3'}\tna\tna\t+\t$dist\t";
				print $pfh "$reg\t$info5[0]\t$regst\t$regen\t+\t$sp\t$DATA{$sp}{'f'}\n" if ($CF{'still_align'});
			} else {					
				if ($gen{$sp}{5}->[3] eq "+") {
					$regst = $gen{$sp}{5}->[1];
					$regen = $gen{$sp}{3}->[2];
					
				} else { #5' will be at the end, 3' at the begining
					$regst = $gen{$sp}{3}->[1];
					$regen = $gen{$sp}{5}->[2];
				}
				$dist = $regen - $regst +1;	
				#OKregions file
				print $ffh "5\t$gen{$sp}{'5'}->[0]\t$gen{$sp}{'5'}->[1]\t$gen{$sp}{'5'}->[2]\t$gen{$sp}{'5'}->[6]\t$gen{$sp}{'5'}->[5]\t";
				print $ffh "3\t$gen{$sp}{'3'}->[0]\t$gen{$sp}{'3'}->[1]\t$gen{$sp}{'3'}->[2]\t$gen{$sp}{'3'}->[6]\t$gen{$sp}{'3'}->[5]\t";
				print $ffh "$gen{$sp}{5}->[3]\t$dist";
				#posi file to extract
				print $pfh "$reg\t$gen{$sp}{'5'}->[0]\t$regst\t$regen\t$gen{$sp}{5}->[3]\t$sp\t$DATA{$sp}{'f'}\n" if ($CF{'still_align'});
			}
			print $fh "$dist\n";
			close $fh;
		}	
		print $ffh "\n";
	}
	close $ffh;
	close $pfh;	
	print $RFH "     => failed checks = $FAILED_CHECKS\n";
	print $RFH "     => OK regions = $ok (tot round $C = $OK)\n" unless ($CF{'file'});
	return 1;			
}

#------------------------------------------------------------------------------
sub check_anchor_hits {
	my $gen = shift;
	# 0		  1		   2      3		 4	  5		 6
	#$t_name,$t_start,$t_end,$strand,$gen,$ratio,$score);
	foreach my $sp (keys %{$gen}) {
		#Check if one anchor doesn't have a hit
		if ((! $gen->{$sp}{5}[0]) || (! $gen->{$sp}{3}[0])) {
			return ("yes","anchor(s) with no hit"); 
		}
		#Check if the target is OK
		if ($gen->{$sp}{5}[0] ne $gen->{$sp}{3}[0]) {
			return ("yes","anchors hit on different targets"); 
		}
		#Check anchor sizes < $CF{'a_max_len'} (150%)
		if (($gen->{$sp}{5}[2] - $gen->{$sp}{5}[1] +1 > $CF{'a_max_len'}) || ($gen->{$sp}{3}[2] - $gen->{$sp}{3}[1] +1 > $CF{'a_max_len'})) {
			return ("yes","anchor size > a_max_len");
		}
		#Check if the score is highest enough
		if ($gen->{$sp}{5}[5] > $CF{'maxratio'} || $gen->{$sp}{3}[5] > $CF{'maxratio'}) {
			return ("yes","second highest hit too close to highest hit")  
		}
		#Check strand, need to be the same (otherwise, inversion)
		if ($gen->{$sp}{5}[3] ne $gen->{$sp}{3}[3]) {
			return ("yes","anchors hit on different strands");
		}
		#Check if weird rearrangement that placed the anchors int he wrong order
		if ($gen->{$sp}{5}[3] eq "+") {
			if ($gen->{$sp}{5}[1] < $gen->{$sp}{5}[2] && $gen->{$sp}{5}[1] > $gen->{$sp}{3}[1]) {
				return ("yes","anchors hit on the plus strand, but they are in the wrong order");
			}		
		} else {
			if ($gen->{$sp}{3}[1] < $gen->{$sp}{3}[2] && $gen->{$sp}{3}[1] > $gen->{$sp}{5}[1]) {
				return ("yes","anchors hit on the minus strand, but they are in the wrong order (possibly at least 2 inversions)");
			}	
		}
		#Check if overlapping anchors
		my $dist;
		if ($gen->{$sp}{5}[3] eq "+") {
			$dist = $gen->{$sp}{3}[1] - $gen->{$sp}{5}[2] + 1;
		} else {
			$dist = $gen->{$sp}{5}[1] - $gen->{$sp}{3}[2] + 1;
		}
		if ($dist < 0) {
			return ("yes","overlaping anchors");
		}
		#Check assembly gaps
		my $chr = $gen->{$sp}{5}[0];
		my $ifnext;
		if ($gen->{$sp}{5}[3] eq "+") {
			$ifnext  = check_overlap_gap($GAPS{$sp}{$chr},$gen->{$sp}{3}[2],$gen->{$sp}{5}[1]);
		} else {
			$ifnext  = check_overlap_gap($GAPS{$sp}{$chr},$gen->{$sp}{5}[1],$gen->{$sp}{3}[2]);
		}
		if ($ifnext eq "yes") {
			return ("yes","overlap with assembly gaps");	
		}
	}
	#if no return yet, region passed filters
	return("no","na");
}


#------------------------------------------------------------------------------
# EXTRACT AND ALIGN
#------------------------------------------------------------------------------
sub extract_regions {
	print $AFH "-------------------------------------------------------------------------------\n";
	print $AFH " ROUND $C\n";
	print $AFH "-------------------------------------------------------------------------------\n";
	print $AFH " --- Extracting regions\n";
	my %reg = ();
	my %noaln = ();
	my %check = ();
	open(my $fh,"<",$POSI) or confess "\t    ERROR - could not open to read $POSI $!\n";
	while(defined(my $l = <$fh>)) {	
		chomp($l);	
		next if ($l !~ /\w/ || $l =~ /^#/);	#avoid blank lines and header lines
		my ($reg,$chr,$st,$en,$strand,$sp,$fa) = split(/\t/,$l);
		unless ($reg && $chr && $st && $en && $strand && $sp && $fa) {
			print $AFH "\t    ERROR - some info are missing ($l)\n";
			next;
		}
		$strand="RC" if ($strand eq "-");
		my $seq = "$ALIGN/$reg.fa";
		
		#get db
		my $db = Bio::DB::Fasta->new($fa) or confess "\t    ERROR - Failed to create Bio::DB::Fasta object from $fa $!\n";
		
		#check existence
		if (-f $seq && ! $check{$reg}) {#ie if file exists despite fact this is the first extraction for this region
			print $AFH "\t    WARN - $seq exists => no extraction, but will still be aligned\n";
			$reg{$seq}=1;
			next; 
		}
		$check{$reg}=1;

		#deal with the | in scaffold names, problem to match
		my $expr = $chr;
		$expr =~ s/\|/\\\|/g;
	
		#name of the sequence extracted will be only gb when there is gi and gb in name (otherwise too long)
		my $chrmod = $chr;
		$chrmod =~ s/^gi\|.*\|gb/gb/;
	
		# prepare new name
		my $newid = ($sp."_".$chrmod.":".$st."-".$en."_".$strand);
	
		#now extract
		my $out = Bio::SeqIO->newFh(-format => 'Fasta', -file=>">>$seq") or confess "\t    ERROR - Failed to create SeqIO FH object from $seq $!\n";
		my $sub = $db->seq($chr,$st,$en) or warn "\t    ERROR - could not extract sequence $chr,$st,$en from $fa\n";
		next if (! $sub);
		my $seqobj = Bio::Seq->new( -display_id => $newid, -seq => $sub) or warn "\t    ERROR - could not create Bio::Seq object for $chr,$st,$en in $fa\n";
		next if (! $seqobj);
		my $len = $seqobj->length;
#		print $AFH "\t      $reg - $sp => $chr:$st-$en => $len\n";
		if ($len > $CF{'r_max_len'}) {
			print $AFH "\t    WARN - $reg.fa won't be aligned (its length = $len nt > max = $CF{'r_max_len'} nt)\n";
			$noaln{$seq}=1;
		} else {
			$seqobj = $seqobj->revcom if ($strand eq "RC"); #reverse complement if hit was antisense
			print $out $seqobj;
			$reg{$seq}=1;
			$KALIGN{$seq}=1 if ($len > $CF{'max_muscle'});
		}	
	}
	close $fh;
	return (\%reg,\%noaln);
}

#------------------------------------------------------------------------------
sub align_regions {
	my $reg = shift;	
	my $noaln = shift;
	print $AFH " --- Aligning regions\n";
	foreach my $seq (keys %{$reg}) {
		print $AFH "       -> $seq\n";
		my $out = "$seq.align.fa";
		if (-f $out || $noaln->{$seq}){
			print $AFH "\t    ERROR - $out exists or had non extracted sequences => Not realigning $seq\n";
		} else {	
			my $alog = "$seq.align.log";
			if ($CF{'alnsoft'} =~ /[Mm]uscle/ && ! $KALIGN{$seq}) {
				print $AFH "          $CF{'alnsoft'} -in $seq -out $out -log $alog -quiet -verbose\n";
				`$CF{'alnsoft'} -in $seq -out $out -log $alog -quiet -verbose 2>> $ALOG`;
				print $AFH "          WARN: muscle failed, Kalign is used instead\n" if (! -f $out);
			}
			if ($KALIGN{$seq} || $CF{'alnsoft'} =~ /[Kk]align/ || ! -f $out) {	
				#if Kalign set for this seq, or in general anyway, or muscle crashed (no out file)	
				print $AFH "          $CF{'alnsoft'} -quiet -i $seq -o $out\n";
				`$CF{'alnsoft'} -quiet -i $seq -o $out 2>> $ALOG`;
				#defaults are: gap open penalty      -gpo   = 6
				#			   gap extension penalty -gpe   = 0.9
			}
			#Easier to look at by eye if re written on 1 line - keep original files for now			
			if (-f $out) {
				unlink $alog unless ($DEBUG);
# 				#perl -pE 'if ($_ !~ /^>/) { s/\n//; } s/>/\n>/;' $out > $seq.align.1l.fa
				`cat $out | perl -pe '/^>/ ? print "\n" : chomp' > $seq.align.1l.fa`;
			} else {
				if ($DEBUG) {
					print $AFH "          WARN: $out missing - check $alog if you want to see why\n";
				} else {
					print $AFH "          WARN: $out missing\n";
				}
			}
		}
	}
	return 1;
}


#------------------------------------------------------------------------------
# AMALYZE GAPS
#------------------------------------------------------------------------------
sub analyze_gaps {
	print $GFH "-------------------------------------------------------------------------------\n";
	print $GFH " ROUND $C\n";
	print $GFH "-------------------------------------------------------------------------------\n";
	print $GFH " --- Analyzing gaps\n";
	
	my @f = `ls $ALIGN/*.align.fa`;
	print $GFH "     Listing gaps for each species\n";
	foreach my $file (@f) {
		chomp($file);
		print $GFH "        -> from $file\n";
		$file =~ s/.*\/(.*)$/$1/;
		extract_and_write_gaps($file);
		#prints small in 1-30 and others in 31-x
	}
	print $GFH "     Concatenating gap files\n";
	concat_gaps("1-30");
	concat_gaps("31-x");
	
	print $GFH "     Obtaining shared gaps\n";
	shared_gaps("1-30",0.80);
	shared_gaps("31-x",0.85);
	
	print $GFH "     Obtaining species specific gaps\n";
	spe_gaps("1-30",0.80);
	spe_gaps("31-x",0.85);
	
	#Then print the data
	print $GFH "     Print the results\n";
	print_gaps("1-30");
	print_gaps("31-x");
	print $GFH " --- Results in: $CURRDEL/__1-30.RESULTS.tab\n";
	print $GFH "                 $CURRDEL/__31-x.RESULTS.tab\n";
	return 1;
}

#----------------------------------------------------------------------------
sub extract_and_write_gaps {
	my $file = shift;	
	my $alignio = Bio::AlignIO->new(-file => $ALIGN."/".$file, -format => 'fasta');
	while(my $aln = $alignio->next_aln()){
		my %seqs = ();
		my $length;
		foreach my $seq ($aln->each_seq() ) {
			my @name = split(/\_/,$seq->display_id);
			my $sp = $name[0];	
			my @sequence = split(//,$seq->seq);
			$seqs{$sp} = \@sequence;
			$length = $seq->length;
			$TOTLEN{$sp}+=$length;
		}
		# Get gaps and print their coordinates - unless the files exist
		my %start = ();
		my %gap_nb = ();
		my %skipped = ();
		for (my $n=1;$n<$length;$n++) {
			foreach my $sp (@SPIDS) {
				#nt for each species at each position is ($seqs{$sp}->[$n]
				#print in a file gap start-end couples, one per line, one file per species
				my $sgaps = $GAPS."/".$file."$sp.gaps.1-30.bed";
				open(my $sfh,">>",$sgaps) or confess "    ERROR - can not open file $sgaps $!";
				my $lgaps = $GAPS."/".$file."$sp.gaps.31-x.bed";
				open(my $lfh,">>",$lgaps) or confess "    ERROR - can not open file $lgaps $!";
				my $reg = $file; #reg4-164.fa.align.fa
				$reg =~ s/\.fa\.align\.fa$//;
				if ($seqs{$sp}->[$n] eq "-" && $seqs{$sp}->[$n-1] ne "-" && ($n > 1 && $seqs{$sp}->[$n-2] ne "-")) { #this is a gap opening
				# "-" at position $n => this is a gap opening; IF $n-1 ne "-" AND IF $n-2 ne "-" [ie situation TGC-----A-----TGC]
					$start{$sp} = $n+1;
					$skipped{$sp}++;
				}
				if ($seqs{$sp}->[$n-1] eq "-" && $seqs{$sp}->[$n] ne "-" && ($start{$sp} && $seqs{$sp}->[$n+1] eq "-")) { 
				#this is a gap ending; unless gap is on begining of alignement - I don't want to count these since I don't know their length
				#and unless it's a case like A in TGC-----A-----TGC (if one or more was seen, start coordinate was "corrected")
					my $end = $n+1;
					my $len = $end - $start{$sp};
					if ($skipped{$sp}) {
						$len = $len - $skipped{$sp}; #correct length
						$skipped{$sp} = 0;
					}	
					if ($len > 30) {
						$gap_nb{$sp}{'l'}++;
						print $lfh "$reg\t$start{$sp}\t$end\t$gap_nb{$sp}{'l'}\t.\t+\t$len\t$length\n";
					} else {
						$gap_nb{$sp}{'s'}++;
						print $sfh "$reg\t$start{$sp}\t$end\t$gap_nb{$sp}{'s'}\t.\t+\t$len\t$length\n";
					}
				}
				close $lfh;
				close $sfh;
			}
		}
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub concat_gaps {
	my $t = shift;
	foreach my $sp (@SPIDS) {
		my $out = "$GAPS/$sp.gaps.$t.bed";
		`cat $GAPS/*$sp.gaps.$t.bed > $out`;
	}
	return 1;
}

#-------------------------------------------------------------------------------
sub shared_gaps {
	my $type = shift;
	my $f = shift;
	#loop through the levels
	foreach my $lvl (sort { $b <=> $a } keys %NWCK) { 
		my @sp = @{$SPIDS{$lvl}};
		my $par = $PAR{$lvl};
		if ($lvl == $MAX) {
			print $GFH "     GROUP $lvl => @sp (root)\n";	
			#top level is the "all species" => all shared gaps => not orientable		
			$ANC = "gaps.$type";
			$CURR = $lvl;
			split_gaps($SPIDS{$lvl},"no",$type,$f); 		
		} else {
			my @parsp = @{$SPIDS{$par}};
			print $GFH "     GROUP $lvl => @sp (parent = $par => @parsp)\n";	
			#Each time, set $ANC and $CURR, with
			#   $ANC = suffix file name from the first round of "split_gaps" of the last common ancestor group processed
			#   $CURR = label of this group in output file names => label the species in the list of interest
			#Use the parent for the "not shared"
			$ANC = "gaps.$type.not-shared.".$par;
			$CURR = $lvl;
			#the use of 'minus' avoids having to retype species names all the time; second list is subtracted from the first one
			my @not_sp = minus($SPIDS{$MAX},$SPIDS{$lvl});		
			split_gaps($SPIDS{$lvl},\@not_sp,$type,$f);
		}
	}
	return 1;
}

#----------------------------------------------------------------------------
sub split_gaps {
	my ($to_loop,$to_sub,$type,$f) = @_;
	my @to_loop = @{$to_loop};
	my $ok = 1;
	my %noshared = (); #check hash for shared gaps, if none
	foreach my $specie (@to_loop) {
		print $GFH "        $specie, gaps: $type nt\n";
		#DEFINE FILE NAMES
		my $shared     = "$GAPS/$specie.gaps.$type.shared.$CURR.bed"; #final output will shared gaps
		my $not_shared = "$GAPS/$specie.gaps.$type.not-shared.$CURR.bed"; #new output with these gaps above subtracted
		#LOOP
		my $sp = 0;
		my $i_in;
		if ($to_sub eq "no") {
			$i_in = "$GAPS/$specie.$ANC.bed"; #all gaps in this case
		} else {
			$i_in = "$GAPS/$specie.$ANC.bed"; #will be longer but now hard to automatize what group is the "parent"
		}
		my $i_in_temp = $i_in; #to avoid rw of $i_in during the looping stuff
		print $GFH "              (shared gaps in $shared; not shared in $not_shared)\n" if (! $noshared{$CURR});
		SP: foreach my $othersp (@to_loop) { 
			my $i_out_temp;
			if ($othersp ne $specie) { #useless to intersect if same species
				$i_out_temp = "$shared.temp.$sp";
				my $i_to_check = "$GAPS/$othersp.$ANC.bed";
				#Avoid doing the intersections if no gap shared between all - but need the files
				if ($noshared{$CURR}) {
					$i_in_temp = $i_out_temp;
					`touch $i_in_temp`;
					last SP;
				} else {
					print $GFH "              intersectBed -a $i_in_temp -b $i_to_check -f $f -r -wa > $i_out_temp\n";
					`$CF{'bedtools'}/intersectBed -a $i_in_temp -b $i_to_check -f $f -r -wa 1> $i_out_temp 2>> $GLOG`;
					# now this output is the input for next round
					$i_in_temp = $i_out_temp;
					#check this i_in_temp file; if empty, it means there won't be shared gaps between the species
					#So no need to intersectBed with the others
					my $lines = 0;
					open(my $prev,"<", $i_in_temp) or warn "           ERROR - can not open to read file $i_in_temp $!";
					while(<$prev>) {
						$lines++;
						last if ($lines > 0);
					}
					close $prev;
					$ok = 0 if ($lines < 1);
					if ($ok == 0) {
						print $GFH "           WARN: No shared gaps between $specie and all species of the group $CURR\n";
						print $GFH "           => skipping the rest of the intersections\n";
						$noshared{$CURR}=1;
						last SP;
					}
				}	
			}
			$sp++;	
		}
		# If there are shared ones:
		# - keep file with shared gaps between all species
		`mv $i_in_temp $shared`; # rename to keep the last file
		# - GET FILE FOR NEXT STEPS => SUBTRACT SHARED STUFF TO GET NON SHARED STUFF (will be input for next time this subroutine is used)
		print $GFH "           -> getting input for the next round\n";
		print $GFH "              subtractBed -a $i_in -b $shared > $not_shared\n";
		`$CF{'bedtools'}/subtractBed -a $i_in -b $shared 1> $not_shared 2>> $GLOG`; #no need flexibility here, same spec => same coords	

		# - NOW FILTER TO GET SHARED THAT ARE SPECIFIC (super stringent since here can't be convergence) unless no need
		#   here I could sub the original files, but would be longer. Only do that if there are no "not-shared" for that species (out of groups)
		unless ($to_sub eq "no" || $noshared{$CURR}) {
			my $i = 0;
			my $i_in_temp2 = $shared; #to avoid rw of $i_in
			print $GFH "           -> loop in complementary list of species to remove gaps shared with these\n";
			foreach my $spec (@{$to_sub}) {
				my $i_to_sub = "$GAPS/$spec.$ANC.bed"; #stuff are going to be subtracted
				$i_to_sub = "$GAPS/$spec.gaps.$type.bed" if (! -e "$GAPS/$spec.$ANC.bed"); #sub the original gaps if needed
				my $i_out_temp2 = "$shared.temp.$i"; #files previously generated = shared but not filtered
				print $GFH "              intersectBed -a $i_in_temp2 -b $i_to_sub -f $f -r -v -wa > $i_out_temp2\n";
				`$CF{'bedtools'}/intersectBed -a $i_in_temp2 -b $i_to_sub -f $f -r -v -wa 1> $i_out_temp2 2>> $GLOG`; 
				     #flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$i_in_temp2 = $i_out_temp2; #now this output is the input for next round
				$i++;
			}
			#      keep file with shared gaps between all species
			`mv $i_in_temp2 $shared` unless ($i_in_temp2 eq $shared); # rename to keep the last file if needed (= if several files)
			#`rm -f $shared.temp*`; # remove temp intermediate files	
		} 
		push(@GAPFILES,$shared);
	}
	`rm -f $GAPS/*.temp*`; # remove temp intermediate files
	return 1;
}

#----------------------------------------------------------------------------
sub spe_gaps {
	my $type = shift;
	my $f = shift;
	foreach my $sp (@SPIDS) {
		print $GFH "        $sp, gaps: $type nt\n";
		my $i_out = "$GAPS/$sp.gaps.$type.specific.bed";
		my $i_in  = "$GAPS/$sp.gaps.$type.bed"; #files previously generated; all gaps for this species
	
		#LOOP to subtract all other gaps
		my $i = 0;
		my $i_in_temp = $i_in; #to avoid rw of $i_in
		foreach my $othersp (@SPIDS) {
			unless ($othersp eq $sp) {
				my $i_to_sub = "$GAPS/$othersp.gaps.$type.bed"; #here I can use all gaps, doesn't change anything
				my $i_out_temp = "$i_out.temp.$i";
				print $GFH "              intersectBed -a $i_in_temp -b $i_to_sub -f $f -r -wa -v > $i_out_temp\n";
				`$CF{'bedtools'}/intersectBed -a $i_in_temp -b $i_to_sub -f $f -r -wa -v 1> $i_out_temp 2>> $GLOG`; 
			     #=> flexibility here, b/c not same species. This syntax means subtracting b from a, return a.
				$i_in_temp = $i_out_temp; #now this output is the input for next round
				$i++;
			}	
		}
		#      keep file with only gaps that are specific to current species
		system "mv $i_in_temp $i_out"; # rename to keep the last file
		system "rm -f $GAPS/$i_out.temp*"; # remove temp intermediate files
		push(@GAPFILES,$i_out);
	}
}

#----------------------------------------------------------------------------
sub print_gaps {
	my $type = shift;
	print $GFH "     Analyzing shared/specific gap files for gaps: $type nt\n";
	my $out = "$CURRDEL/__$type.RESULTS.tab";
	open(my $fh,">$out") or confess "     ERROR - can not open to write $out $!";
	print $fh "#Results from $SCRIPTNAME v$VERSION - ran with $CONFIG_FILE\n";
	print $fh "# => ROUND $C\n";
	print $fh "\n";
	print $fh "#Species groups are as follow:\n";
	print $fh "#.\tid\tspecies_list\n";
	foreach my $lvl (sort { $b <=> $a } keys %NWCK) { 
		my @sp = @{$SPIDS{$lvl}};
		print $fh "group:\t$lvl\t@sp\n";	
	}
	print $fh "\n";
	print $fh "#file_name\tspecies\ttype\tgroup\tcount\tamount\taln_len\t%_of_aln\tdel_nt/10kb\tif_ougroup\n";
	foreach my $file (@GAPFILES) {
		next unless ($file =~ /$type/);
		print $GFH "        -> $file\n";
		#0				1		2	3	4	5	6			7
		#chr1_block.4	start	end	nb	.	+	gap_len		aln_len(block)
		my $fname = $file;
		$fname =~ s/.*\/(.*)$/$1/;
		my @n = split(/\./,$fname); #spid.gaps.1-30.shared.group.bed
		$n[4] = "." if ($n[4] eq "bed");
		if (-f $file) {
			open(my $ifh, "<",$file) or confess "     ERROR - can not open to read $file $!";
			my $totgaplen = 0;
			my %gapcount;
			while (<$ifh>) {
				chomp (my $line = $_);
				my @line = split(/\t/,$line);
				my $gaplen = $line[2]-$line[1];
				my $id = $line[0]."_".$line[3];
				$totgaplen += $gaplen;
				$gapcount{$id}++;
			}
			close $ifh;
			my $gaps = 0;
			foreach my $gap (sort keys %gapcount)  {
				$gaps++;
			}
			my $gap_per = 0;
			$gap_per = $totgaplen / $TOTLEN{$n[0]} * 100;
			my $del_nt = $totgaplen / $TOTLEN{$n[0]} * 10000;
			print $fh "$fname\t$n[0]\t$n[3]\t$n[4]\t$gaps\t$totgaplen\t$TOTLEN{$n[0]}\t$gap_per\t$del_nt";	
		} else {
			print $fh "$fname\t$n[0]\t$n[3]\t$n[4]\t0\t0\t$TOTLEN{$n[0]}\t0\t0";
		}
		#Now add a warning if outgroup
		if ($n[0] eq $OG || ($n[4] ne "." && $n[4] == $MAX)) {
			print $fh "\tOUTGROUP:can't_use_the_numbers\n";	
		} else {
			print $fh "\t.\n";	
		}
	}	
	close $fh;
	return;
}


	
