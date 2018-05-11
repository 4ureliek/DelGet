# DelGet

### SYNOPSIS: 
This code will find micro- and midsize deletions (using gaps of 1-30 nt and >30 nt respectively) between orthologous regions of 3 or more species.

The two final outputs (for 1-30 nt and >30 nt) will contain these columns:

	count	amount	aln_len	%_of_aln	del_nt/10kb
	
A rate can be obtained for the branch by dividing 'del_nt/10kb' by the corresponding My.

#### Important:	

See the v4.6 directory for the code and readme of the version used to generate the data in [Kapusta, Suh and Feschotte (2017) PNAS](http://www.pnas.org/content/114/8/E1460.full) ([see also on BioRxiv](http://biorxiv.org/content/early/2016/10/16/081307))
	
This current version now supports more than 3 species, which means that the gap analysis has changed (it is the same as the one used for [MAFmicrodel](https://github.com/4ureliek/MAF_parsing/tree/master/MAFmicrodel))

Note that the new gap analysis affects cases like this:

	outgroup ACCGTGTGTATGTGTGTGTGCGTGCGCGCGTATGTGTCTGTCTGTGTGCGTGTCTGTACGTGTATATAT
	species1 ACTGTGTGTGTGTGTGTGT------------------------------GTGTCTGTGCGTGTATATAT
	species2 ACTGTGTGTATGTGTGTGTG------------TGTGTGTGTGTGTCTGTGTGTCTGTGCGTGTATATAT
     
     With v4.6, species2 has no species-specific gap and species1 has two (1nt and 17nt); 
        the shared portion is ignored (considered as 'shared').
     With v5+, both species 'keep' the full length gaps as specific, which is more 'biological'.
     If this is a real biological gap, v4.6 would underestimate deletion rates; after comparing 
        same trios of species with v4.6 and v5.3, it looks like this mostly affects microdeletion rates.
     This could be a technical issue (sequencing of a simple repeat), but since these regions are 
        also more prone to indels, these could be real.
     Thus, it is good to keep this type of examples in mind while discussing microdeletion rates.  
	
### USAGE: 
 Regions analyzed can be:
   1) loaded from a file (only Repeat Masker .out files for now) 
   2) chosen randomly in one of the species (see config files for numbers, length, etc)
   
 All options and paths are set in a config file, that can be obtained by running: 
 
	perl DelGet.pl -c configfile.txt
    
 This is were you NEED to define parameters, software and file locations, etc, by editing any line that does not start with a # and has a =.
      
 For more info about the config file and the outputs and/or if you have never used this pipeline, type:

	perl DelGet.pl -h

 Once the config file is edited, typically:
 
	 nohup perl $scriptname DelGet.manipname.txt &

 Sometimes the code stalls - if no 'recent' Del/Deletion.X folders, and if the only job still running is $SCRIPTNAME DelGet.manipname.txt (or if it has been running for too long), consider killing the talled job and restarting, with:
 
     nohup perl $SCRIPTNAME --restart DelGet.manipname.txt &
	
### REQUIREMENTS:
	(already in Lib) Perl Array::Unique = https://metacpan.org/source/SZABGAB/Array-Unique-0.08/lib/Array/Unique.pm
	(already in Lib) Perl Array::Transpose::Ragged = https://metacpan.org/pod/Array::Transpose::Ragged 
	Blat = http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3   
	Muscle = http://www.drive5.com/muscle/downloads.htm   
	Kalign = http://msa.sbc.su.se/cgi-bin/msa.cgi?mode=downloads (if large regions to align, muscle crashes)

### CITATION:
[Kapusta, Suh and Feschotte (2017) PNAS (doi: 10.1073/pnas.1616702114)](http://www.pnas.org/content/114/8/E1460.full) and DelGet.pl vX.X, available at https://github.com/4ureliek/DelGet

### OUTPUTS
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


### STEPS: 
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

