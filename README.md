# DelGet
Get micro and medium deletions specific to 2 closely related species after their split.
All options and paths are set in the file DelGet--0--CONFIG.pl (simple text file)

Usage     : perl DelGet.pl <config_file>
Typically : perl DelGet.pl DelGet--0--CONFIG_manipname.pl

	
PURPOSE : 
   General question = get medium size deletion rates
   Requires the UCSC BLAT stand alone software, see http://genome.ucsc.edu/FAQ/FAQblat.htmlblat3
   Written for 3 species (referred as 1, 2 and 3)
      1) From a \"reference\" genome (=1), get random regions as: [anchor1]=====Genome1.region=====[anchor2]
      2) Blat anchor1 and anchor2 sequences agains genome2 and genome3
      3) Get length between anchors: [anchor1]=====Genome2.region=====[anchor2] and [anchor1]=====Genome3.region=====[anchor2]

   NOTE THAT:
      - anchors are chosen based on their score + if second highest score is < XX% of the highest score (see \"maxratio\" variable in CONFIG file)
      - anchors need to be on same target sequence in both genomes
      - anchors on genome2 and genome3 need to be < XX nt ("a_max_len" variable in CONFIG file);
      - anchor1 and anchor2 hits need to be on same strand
      - and obviously no assembly gaps in any of the regions
				
CONFIG file :
   This is were you NEED to define the 3 genome locations + gap files, by editing path between the quotes.
      These gap files list assembly gaps coordinates. Files need to be as UCSC format, and obvisouly name of sequences need to match names of sequences in your genome file. 
      If no gap file is available on UCSC for your assembly, use the provided fasta_get_gaps.pl
      If you never see in the log file \"gap in this region in gen1 => next random posi\" or \"there are gaps in this region\", then there is an issue with these files, please check format etc.
   This is were you NEED to define BLAT software location
   This is were you NEED to define the pipeline parameters
   Because of 2 different indexing happen, you need to define different location for genomes to extract sequences from. Sym links work.
		

OUTPUTS :   
   0) From this pipeline script
      -> _DelGet.log
         You WANT to check it. Steps are detailed, as well as files names, checking steps etc.
	
   1) FROM SCRIPT DelGet--1--get_regions.pl
      INSIDE FOLDER Deletions.$c ($c being a counter to avoid erasing previous results since script will be run several times)
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
         Contains specific gaps. See script itself for more details.
