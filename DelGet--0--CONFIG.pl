####################################################################################################################
# CONFIGURATION FILE FOR: DelGet.pl
# version of the pipeline = 3.7
#
# if you never used this script, type: perl DelGet.pl -h 
#
####################################################################################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
####################################################################################################################
# please edit paths or values. This text file is saved as a .pl only in order to get nice colors when opened in Notepad++.
# do not remove the "#" signs, as they mean the text shouldn't be interpreted (comments) by perl
# do not put "/" at the end of the paths
# keep the quotes around the paths
# the ";" at the end is mandatory


####################################################################################################################
# VARIABLES / PARAMETERS OF THE RUN
####################################################################################################################
# FOR SCRIPT --1--
# NUMBER OF TOTAL REGIONS = to stop the loop
# 1) number of random positions that will be treated per round (sequence extracted, blat against gen2 and gen3, checking steps etc).
  $a_per_round = 500; #[blat step is a limitation in the script and requires to put genome in memory every time - don't put that number too low or too high]
# 2) TOTAL number_of_randomization 
  $randtot = 100000; #Loop inside the pipeline; rounds will be done $randnb by $randnb anyway
# 3) Number per run that need to be successful
  $randnb = 100; #[this is arbitrarily chosen based on when script was stopping to write outputs for some groups of species]
  
# length beetween anchors  
  $anch_dist = 25000; #=distance in genome1 between the 2 anchors [anchor1]<----XXnt---->[anchor2]

# length of anchors.
  $anch_len = 100; #100 is a good number.
  
# Highest ratio (highest_score / second_highest_score). This will allow to filter for hits "not high enough"
  $maxratio = 0.9; #let's say we need highest score / second score < 0.9 (ie second score is max 90% of highest)
  
# Minimum length of a N strech for it to be considered as an assembly gap
  $mingaplen = 50; #[typically assembly gaps are 100nt]

# total length of regions (+2nt) [ADAPT ONLY IF WEIRD ANCHOR HITS]
  $a_max_len = $anch_len + 50*($anch_len)/100; #[this is to avoid huge span of the hit of anchors in other genome]
  
# [DO NOT CHANGE] total length of regions (+2nt)
  $minlen = $anch_dist + 2*$anch_len;
  
  
# FOR SCRIPT --2--get-sev-seq_align
 # this is to set max length of sequence to extract to align, avoid crashing muscle
  $multip = 7 if ($anch_dist <= 10000); #10kb intially => 70kb max [muscle].
  $multip = 4 if (($anch_dist <= 50000) && ($anch_dist > 10000)); #25kb => 100kb max [muscle], 50kb => 200kb max [kalign software].


####################################################################################################################
# PATHS
####################################################################################################################
# FOR PIPELINE / SEVERAL SCRIPTS
##########################################################
# path to other scripts required
  $script_path = "/home/aurelie/bin/my_Deletions"; #erv
  
# path of the folder that will contain all outputfiles. "." means where pipeline is started [a folder will be created]
  $path = ".";
  
# genome IDs
# these IDs will be used to name the files and added in extracted sequences names. Don't chose something too long. 4 letters of species IDs are the best.
  $IDgen1 = "Rmac";
  $IDgen2 = "Hsap";
  $IDgen3 = "Ptro";
  
##########################################################
# FOR SCRIPT --1--
##########################################################
# blat software
  $BLATSOFT = "/home/software/ucsc/blat/blat"; #erv

# path to a file with OK regions ("_OKregions.tab", output of a previous run), to be loaded and avoid overlap 
# Define path of previous output (next outputs will be automatically concatenated to this one), needed only if not in Deletions.X folder.
#	=> IF PREVIOUS OUTPUTS ARE IN FOLDERS Deletions.X, DO NOT SPECIFY PATH, IT WILL BE ADDED
#	   This means that if you don't want to load regions of previous runs, change names or move previous folders
  # $OKregions = "/data/aurelie/GenomeSize/Deletions/XXX/_OKregions.tab";
  $OKregions = "no"; 
  $if_OKregions = "yes"; # => if OKregions is = "no", path of previous output will be detected automatically and all previous regions will be concatenated in a file if relevant
					     #    this file will be created in the script and be the result of: concatenation of Run.(-)1/_OKregions.tab and Run.(-)2/_OKregions.all.tab [if relevant]
					     # => if OKregions is = "/path_to/file.tab", then it will be used as the first previous output and only at second run it will be automatic 
  # $if_OKregions = "no"; # means that no previous output should be loaded AT ANY TIME [not recommended]

# genomes files
# note that you need writing access over there, to create index files
  #erv
  $genone   =  "/data/genomes/Rhesus/rheMac3.fa";
  $gentwo   =  "/data/genomes/human/hg19/chr/hg19.chr.fa";
  $genthree =  "/data/genomes/Chimpanzee/panTro4/panTro4.fa";


# gap coordinates files
# These gap files list assembly gaps coordinates. Files need to be as UCSC format, and obvisouly name of sequences need to match names of sequences in your genome file. 
# If no gap file available on UCSC for your assembly, use the provided script fasta_get_gaps.pl 
  #erv
  $gapone   =  "/data/aurelie/GenomeSize/Deletions/Primates/PrimatesGaps/rheMac3.gaps.tab";
  $gaptwo   =  "/data/aurelie/GenomeSize/Deletions/Primates/PrimatesGaps/hg19.gaps.tab";
  $gapthree =  "/data/aurelie/GenomeSize/Deletions/Primates/PrimatesGaps/panTro4.gaps.tab";

  
##########################################################
# FOR SCRIPT --2--get-sev-seq_align
##########################################################
# Kalign software [needed only for large alignments], see http://www.biomedcentral.com/1471-2105/6/298/
  $KALIGNSOFT = "/home/software/Kalign2/kalign"; #erv
  
# muscle software  
  $MUSCLESOFT = "/home/software/muscle3.8.31/muscle3.8.31"; #erv

# genome files TO EXTRACT SEQUENCES = needed to write in posifile
  #erv
  $genone_a   =  "/data/aurelie/GenomeSize/Deletions/genomes_toalign/rheMac3.fa";
  $gentwo_a   =  "/data/aurelie/GenomeSize/Deletions/genomes_toalign/hg19.chr.fa";
  $genthree_a =  "/data/aurelie/GenomeSize/Deletions/genomes_toalign/panTro4.fa";
  
# [DO NOT CHANGE]: decide aln software
  ($anch_dist >= 40000)?($ALNSOFT = $KALIGNSOFT):($ALNSOFT = $MUSCLESOFT);
  
##########################################################
# FOR MASKING
##########################################################
# repeat masker software
  $RMSOFT = "/home/software/RepeatMasker/RepeatMasker"; #erv

# masking librarie(s)
# If blank, no masking will be done. Several libraries can be defined and a folder will be created using the @masking_folders list
# Note that for now, a -lib is needed => use RM script to convert their embl to fasta if all stuff should be used
# erv
  @mask = ("/data/aurelie/GenomeSize/Deletions/Primates/TEdb_to_mask/20130520_repbase_ancientTEs.fa","/data/aurelie/GenomeSize/Deletions/Primates/TEdb_to_mask/RepeatMaskerLib.20130422.fa");
  @masking_folders = ("masked.ancient","masked.RM20130422");
#

#ensure last returned value is true (to load as require in scripts) => DO NOT REMOVE
1;
