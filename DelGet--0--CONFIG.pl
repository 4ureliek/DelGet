####################################################################################################################
# CONFIGURATION FILE FOR: DelGet.pl
# version of the pipeline = 4.2
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
# 1) TOTAL number of randomization 
  $randtot = 10000; 
     # This number will be reached using a loop inside the pipeline, running all scripts
     # Rounds will be done $randnb by $randnb.
# 2) Number per run that need to be successful before going to extracting and aligning sequences
  $randnb = 100;
      # After a while, it is possible that very little or no anchors manage to get through the filters
      # Running in loop allows to still get outputs and results, even if script has to be killed for that reason or any other
      # Use the script in utilities to cat the outputs
# 3) number of random positions that will be treated per round (anchor sequences extracted from outgroup, blat against gen2 and gen3, checking steps etc).
  $a_per_round = 500; 
     # This step is repeated inside script --1-- until $randnb is reached. 
     # Blat step is a limitation in the script and requires to put genome in memory every time;
     # therefore, if species are closely related or if assembly is good, you can lower this number, 
     # but given the usual success rate even in primates, at least $randnb x2 is advised.
     # Do x5 if less good assembly or more distant species, or if they have a lot of recent TEs for example.
  
# length beetween anchors  
  $anch_dist = 25000; #=distance in genome1 between the 2 anchors [anchor1]<----XXnt---->[anchor2]

# length of anchors.
  $anch_len = 100; #100 is a good number.
  
# Highest ratio (highest_score / second_highest_score). This will allow to filter for hits "not high enough"
  $maxratio = 0.9; #let's say we need highest score / second score < 0.9 (ie second score is max 90% of highest)
  
# Minimum length of a N strech for it to be considered as an assembly gap
  $mingaplen = 50; # typically assembly gaps are 100nt 

# total length of regions (+2nt) [ADAPT ONLY IF WEIRD ANCHOR HITS]
  $a_max_len = $anch_len + 50*($anch_len)/100; # This is to avoid huge span of the hit of anchors in other genome
  
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
# path of the folder that will contain all outputfiles. "." means where pipeline is started [a folder will be created]
  $path = "./Del";
  
# genome IDs
# these IDs will be used to name the files and added in extracted sequences names. Don't chose something too long. 4 letters of species IDs are the best.
  $IDgen1 = "Rmac";
  $IDgen2 = "Hsap";
  $IDgen3 = "Ptro";
  
##########################################################
# FOR SCRIPT --1--
##########################################################
# blat software
  $BLATSOFT = "/home/software/ucsc/blat/blat";

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
  $genone   =  "/data/genomes/Rhesus/rheMac3.fa";
  $gentwo   =  "/data/genomes/human/hg19/chr/hg19.fa";
  $genthree =  "/data/genomes/Chimpanzee/panTro4/panTro4.fa";


# gap coordinates files
# These gap files list assembly gaps coordinates. Files need to be as UCSC format, and obvisouly name of sequences need to match names of sequences in your genome file. 
# If no gap file available on UCSC for your assembly, use the provided script fasta_get_gaps.pl 
  $gapone   =  "/data/Deletions/Primates/PrimatesGaps/rheMac3.gaps.tab";
  $gaptwo   =  "/data/Deletions/Primates/PrimatesGaps/hg19.gaps.tab";
  $gapthree =  "/data/Deletions/Primates/PrimatesGaps/panTro4.gaps.tab";

  
##########################################################
# FOR SCRIPT --2--get-sev-seq_align
##########################################################
# Kalign software [needed only for large alignments], see http://www.biomedcentral.com/1471-2105/6/298/
  $kalign = "/home/software/Kalign2/kalign"; #erv
  
# muscle software  
  $muscle = "/home/software/muscle3.8.31/muscle3.8.31"; #erv

# genome files TO EXTRACT SEQUENCES = needed to write in posifile
  #Note that you can simply create symbolic links of the genome files, with ln -s /data/genomes/Rhesus/rheMac3.fa /data/Deletions/genomes_toalign/rheMac3.fa
  $genone_a   =  "/data/Deletions/genomes_toalign/rheMac3.fa";
  $gentwo_a   =  "/data/Deletions/genomes_toalign/hg19.fa";
  $genthree_a =  "/data/Deletions/genomes_toalign/panTro4.fa";
  
# [DO NOT CHANGE]: decide aln software
  ($anch_dist >= 40000)?($ALNSOFT = $kalign):($ALNSOFT = $muscle);
  
##########################################################
# FOR MASKING
##########################################################
# repeat masker software
  $RMSOFT = "/home/software/RepeatMasker/RepeatMasker";

# masking librarie(s)
# If blank, no masking will be done. 
# Several libraries can be defined, just separate by a , in the @mask and the @masking_folders list
# As many masking_folders elements as mask files are required
# -lib is needed in running repeat masker => use Repeat Masker script in their utilities to convert their embl to fasta if you would like to mask with the whole library
  @mask = ("/data/Deletions/Primates/TEdb_to_mask/RepeatMaskerLib.20140131.fa","/data/Deletions/Primates/TEdb_to_mask/TE.subset.fa");
  @masking_folders = ("masked.RM20130422","masked.subset");
#

#ensure last returned value is true (to load as require in scripts) => DO NOT REMOVE
1;
