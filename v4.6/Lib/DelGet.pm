#!/usr/bin/perl -w
#######################################################
# SUBROUTINES
# FOR SUBSET OF SCRIPTS = DelGet.pl
# => defined as DelGet package
######################################################
package DelGet;
my $version = "4.6";

######################################################
# QUITE GENERAL FILE ARCHITECTURE STUFF
######################################################
#----------------------------------------------------------------------------
# get a filename from a full path
# my $genone_name = DelGet::filename($genone);
#----------------------------------------------------------------------------
sub filename {
	my($name) = shift;
	$name =~ s/.*\/(.*)$/$1/;
	return $name;
}

#----------------------------------------------------------------------------
# make output directory
#----------------------------------------------------------------------------
sub make_out_dir {
	my($path) = shift;
	my $c=0;
	until (!-d "$path/Deletions.$c"){
		$c++;
	}
	my $pathtemp = "$path/Deletions.$c";
	mkdir ($pathtemp, 0755) or die "\t    ERROR: SUB make_out_dir: can not mkdir $pathtemp $!";
	return ($pathtemp,$c);
}

######################################################
# QUITE GENERAL FASTA STUFF
######################################################
#----------------------------------------------------------------------------
# calculate length of the genome, for sequences > $minlen - includes check steps to avoid repeating this if length already calculated
# => my ($totlength) = DelGet::get_tot_len_filtered($genome,$log,$minlen);
# 	 note that log file need to be closed in the main before calling the subroutine, and reopened after
#----------------------------------------------------------------------------
sub get_tot_len_filtered {
	my $genome = shift;
	my $log = shift;
	my $minlen = shift;
	
	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open to write $log $!\n";	

	#calculating total length of the genome (of the database)
	my $lengthfile = "$genome.min-$minlen.length.txt";	
	unless (-e $lengthfile) {
		print $log_fh "\t\tTotal length of genome for scaffolds > $minlen not known ($lengthfile does not exists) - Calculating length...\n";
		
		# index the genome and connect to the fasta file
		my $reindex;
		my $indexfile = "$genome.index";
		if (-e $indexfile) {
			$reindex = 0;
			print $log_fh "\t\tGenome previously indexed ($indexfile exists) - Skipping indexing step...\n";
		} else {
			$reindex = 1;
			print $log_fh "\t\tGenome not indexed ($indexfile does not exists) - Indexing...\n";
		}
		my $db = Bio::DB::Fasta->new($genome, -reindex=>$reindex) or die print $log_fh "\t\tERROR: could not create Bio::DB::Fasta object from $genome $!\n";
		
		#Now loop and get total length
		my @dbIDs = $db->get_all_ids();			
		open (my $len_fh, ">", $lengthfile) or die print $log_fh "\t\t ERROR: could not create file $lengthfile $!\n\n";
		$GenLen = 0;
		foreach my $ID (@dbIDs) {
			my $obj = $db->get_Seq_by_id($ID);
			my $len = $obj->length;
			$GenLen += $len if ($len > ($minlen-2));;
		}
		print $len_fh $GenLen;
		close $len_fh;
	} else {
		print $log_fh "\t\tTotal length of genome for scaffolds > $minlen has been previously calculated ($lengthfile exists)\n";
		open (my $len_fh2, "<", $lengthfile) or die print $log_fh "ERROR: could not open $lengthfile $!\n";
		while (<$len_fh2>) {
			$GenLen = $_;
		}
		close $len_fh2;
	}
	print $log_fh "\t\t => total len = $GenLen nt\n";
	close $log_fh;
	return($GenLen);
}

#----------------------------------------------------------------------------
# Get all lengths of the genome - includes check steps to avoid repeating this if length already calculated
# => my ($geninfos_ref,$dbIDs) = DelGet::get_all_len_filtered($genone,$log,$minlen,$totlength,$infofile,$randtot);
# 	 note that log file need to be closed in the main before calling the subroutine, and reopened after
#----------------------------------------------------------------------------
sub get_all_len_filtered {
	my ($genome,$log,$minlen,$GenLen,$lengthfile,$randtot) = @_;	
	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open to write $log $!\n";		
	# Get lengths and ratios
	print $log_fh "\t--- Getting all lengths and number of random positions per scaffold\n";
	my @dbIDs_filtered = ();
	my $geninfos = ();
	unless (-e $lengthfile) {
		print $log_fh "\t\tAll lengths of scaffolds > $minlen not known ($lengthfile does not exists) => Calculating lengths and getting info in hash\n";
		# index the genome and connect to the fasta file
		my $reindex;
		my $indexfile = "$genome.index";
		if (-e $indexfile) {
			$reindex = 0;
			print $log_fh "\t\tGenome previously indexed ($indexfile exists) - Skipping indexing step...\n";
		} else {
			$reindex = 1;
			print $log_fh "\t\tGenome not indexed ($indexfile does not exists) - Indexing...\n";
		}
		my $db = Bio::DB::Fasta->new($genome, -reindex=>$reindex) or die print $log_fh "\t\tERROR: could not create Bio::DB::Fasta object from $genome $!\n";
		
		#create list of the ID of the genome file + initialize for filtered one
		my @dbIDs = $db->get_all_ids();
		
		open (my $len_fh, ">", $lengthfile) or die print $log_fh "\t\tERROR: could not create file $lengthfile $!\n\n";
		foreach my $ID (@dbIDs) {
			my $obj = $db->get_Seq_by_id($ID);
			my $len = $obj->length;		
			#filter
			if ($len > ($minlen-2)) {
				my $ratio = int($randtot / $GenLen * $len);
				$ratio = 1 if ($ratio == 0); #bug fix 2015 05 13
				print $len_fh "$ID\t$len\t$ratio\n";
				my @infos = ($len,$ratio);
				$geninfos{$ID} = \@infos;
				push(@dbIDs_filtered,$ID);
			}
		}
		close $len_fh;
	} else {
		print $log_fh "\t\tAll lengths of scaffolds > $minlen has been previously calculated ($lengthfile exists) => loading in hash\n";
		open (my $len_fh2, "<", $lengthfile) or die print $log_fh "ERROR: could not open $lengthfile $!\n";
		while (<$len_fh2>) {
			chomp (my $line = $_);
			my ($ID,$len,$ratio) = split (/\t/,$line);
			my @infos = ($len,$ratio);
			$geninfos{$ID} = \@infos;
			push(@dbIDs_filtered,$ID);
		}
		close $len_fh2;
	}
	close $log_fh;	
	return(\%geninfos,\@dbIDs_filtered);
}



######################################################
# QUITE SPECIFIC TO THIS PIPELINE
######################################################
#----------------------------------------------------------------------------
# print config file if asked
# print_config($conf) if ($conf);
#----------------------------------------------------------------------------
sub print_config {
	my $cfile = shift;
	open(my $fh, ">", $cfile) or die "\t    ERROR: SUB print_config: could not open to write $cfile!\n";
	print $fh "####################################################################################################################\n";
	print $fh "# CONFIGURATION FILE FOR: DelGet.pl\n";
	print $fh "# version of the pipeline = $version\n";
	print $fh "# if you never used this script, type: perl DelGet.pl -h \n";
	print $fh "####################################################################################################################\n";
	print $fh "# Author  :  Aurelie Kapusta\n";
	print $fh "# email   :  4urelie.k\@gmail.com\n";
	print $fh "# GitHub  :  https://github.com/4ureliek\n";
	print $fh "####################################################################################################################\n";
	print $fh "# please edit paths or values anytime you see a non commented line and a \"=\"\n";
	print $fh "# do not remove the # signs and do not put / at the end of the paths\n";
	print $fh "\n";
	print $fh "####################################################################################################################\n";
	print $fh "# VARIABLES / PARAMETERS OF THE RUN\n";
	print $fh "####################################################################################################################\n";
	print $fh "# FOR SCRIPT --1--\n";
	print $fh "# NUMBER OF TOTAL REGIONS = to stop the loop\n";
	print $fh "# 1) TOTAL number of randomization \n";
	print $fh "  randtot = 10000\n";
	print $fh "     # This number will be reached using a loop inside the pipeline, running all scripts\n";
	print $fh "     # Rounds will be done randnb by randnb and printed as it progresses\n";
	print $fh "# 2) Number per run that need to be successful before going to extracting and aligning sequences\n";
	print $fh "  randnb = 100\n";
	print $fh "      # After a while, it is possible that very little or no anchors manage to get through the filters\n";
	print $fh "      # Running in loop allows to still get outputs and results, even if script has to be killed for that reason or any other\n";
	print $fh "      # Use the script in utilities to gather the gap lengths from the outputs\n";
	print $fh "# 3) number of random positions that will be treated per round (anchor sequences extracted from outgroup, blat against gen2 and gen3, checking steps etc).\n";
	print $fh "  a_per_round = 500\n";
	print $fh "     # This step is repeated inside script --1-- until randnb is reached. \n";
	print $fh "     # Blat step is a limitation in the script and requires to put genome in memory every time\n";
	print $fh "     # therefore, if species are closely related or if assembly is good, you can lower this number, \n";
	print $fh "     # but given the usual success rate even in primates, at least randnb x2 is advised.\n";
	print $fh "     # Do x5 if less good assembly or more distant species, or if they have a lot of recent TEs for example.\n";
	print $fh "  \n";
	print $fh "# length beetween anchors, in nt\n";
	print $fh "  anch_dist = 25000 #=distance in genome1 between the 2 anchors [anchor1]<----XXnt---->[anchor2]\n";
	print $fh "\n";
	print $fh "# length of anchors, in nt\n";
	print $fh "  anch_len = 100 #100 is a good number; 80 also gave good results\n";
	print $fh "  \n";
	print $fh "# Highest ratio (highest_score / second_highest_score). This will allow to filter for hits not high enough\n";
	print $fh "  maxratio = 0.9 #let's say we need highest score / second score < 0.9 (ie second score is max 90% of highest)\n";
	print $fh "  \n";
	print $fh "# Minimum length of a N strech for it to be considered as an assembly gap\n";
	print $fh "  mingaplen = 50 # typically assembly gaps are 100nt\n";
	print $fh "\n";
	print $fh "# [FYI] this is how the maximum length of the anchor is set (+2nt)\n";
	print $fh " # a_max_len = anch_len + 50*(anch_len)/100 # This is to avoid huge span of the hit of anchors in other genome\n";
	print $fh "  \n";
	print $fh "# [FYI] this is how the total length of regions (+2nt) is calculated\n";
	print $fh " # minlen = anch_dist + 2*anch_len\n";
	print $fh "  \n";
	print $fh "  \n";
	print $fh "# FOR SCRIPT --2--get-sev-seq_align\n";
	print $fh " # [FYI] this how the max length of sequence to extract to align, to avoid crashing muscle\n";
	print $fh " # multip = 7 if anch_dist <= 10000 #10kb intially => 70kb max [muscle]\n";
	print $fh " # multip = 4 if anch_dist <= 50000 and anch_dist > 10000 #25kb => 100kb max [muscle], 50kb => 200kb max [kalign software]\n";
	print $fh "\n";
	print $fh "\n";
	print $fh "####################################################################################################################\n";
	print $fh "# PATHS\n";
	print $fh "####################################################################################################################\n";
	print $fh "# FOR PIPELINE / SEVERAL SCRIPTS\n";
	print $fh "##########################################################\n";
	print $fh "# path of the folder that will contain all outputfiles. \".\" means where pipeline is started [the directory Del will be created]\n";
	print $fh "  path = ./Del\n";
	print $fh "  \n";
	print $fh "# genome IDs\n";
	print $fh "# these IDs will be used to name the files and added in extracted sequences names. \n";
	print $fh "# Don't chose something too long. 4 letters of species IDs are the best, or assembly IDs such as hg38, mm10 etc. \n";
	print $fh "# Species1 is the outgroup\n";
	print $fh "  IDgen1 = Species1 #-------|_\n";
	print $fh "  IDgen2 = Species2 #---|___|\n";
	print $fh "  IDgen3 = Species3 #---|\n";
	print $fh "  \n";
	print $fh "##########################################################\n";
	print $fh "# FOR SCRIPT --1--\n";
	print $fh "##########################################################\n";
	print $fh "# blat software\n";
	print $fh "  BLATSOFT = /home/software/ucsc/blat/blat\n";
	print $fh "\n";
	print $fh "# Behavior for previous outputs (allows to append new regions to existing previous runs):\n";
	print $fh "  OKregions = no\n";
	print $fh "#	   This means that you do not want to load a specific file; but previous regions that are\n";
	print $fh "#	   in the folders Deletions.X in the path defined above will be loaded, to avoid overlaps\n";
	print $fh "#	   This means that if you don't want to load regions of previous runs, change the names or move previous the output folders\n";
	print $fh "#       If you wish to load a file from a different folder, set this OKregions with path to a file with OK regions, such as: OKregions = /data/DelGet/MyRun/PrevDel/_OKregions.tab\n";
	print $fh "   \n";
	print $fh "# Behavior for loading outputs generated during the run (to avoid overlaps):\n";	
	print $fh "  if_OKregions = yes # => if OKregions above is set as no, ALL previous OK regions will be concatenated by the pipeline,\n";
	print $fh "					    #    by concatenating Deletions.Run(-)1/_OKregions.tab and the Deletions.Run(-)2/_OKregions.all.tab (if any)\n";
	print $fh "					    # => if OKregions above is set as a path to a file, then it will be used as the first previous output\n";
	print $fh "                     # Setting this to no (if_OKregions = no) is NOT RECOMMANDED, because it means that no previous output should be loaded AT ANY TIME\n";
	print $fh "\n";
	print $fh "# genomes files (can be symbolic links created with ln -s)\n";
	print $fh "# note that you need writing access over there, to create index files\n";
	print $fh "  genone   =  /data/DelGet/MyRun/data/genome1.fa\n";
	print $fh "  gentwo   =  /data/DelGet/MyRun/data/genome2.fa\n";
	print $fh "  genthree =  /data/DelGet/MyRun/data/genome3.fa\n";
	print $fh "\n";
	print $fh "\n";
	print $fh "# gap coordinates files\n";
	print $fh "# These gap files list assembly gaps coordinates. Files need to be as UCSC format, and obvisouly name of sequences need to match names of sequences in your genome file. \n";
	print $fh "# If no gap file available on UCSC for your assembly, use the provided script fasta_get_gaps.pl \n";
	print $fh "  gapone   =  /data/DelGet/MyRun/data/genome1.gaps.tab\n";
	print $fh "  gaptwo   =  /data/DelGet/MyRun/data/genome2.gaps.tab\n";
	print $fh "  gapthree =  /data/DelGet/MyRun/data/genome3.gaps.tab\n";
	print $fh "\n";
	print $fh "  \n";
	print $fh "##########################################################\n";
	print $fh "# FOR SCRIPT --2--get-sev-seq_align\n";
	print $fh "##########################################################\n";
	print $fh "# Kalign software [needed only for large alignments], see http://www.biomedcentral.com/1471-2105/6/298/\n";
	print $fh "  kalign = /home/software/Kalign2/kalign\n";
	print $fh "  \n";
	print $fh "# muscle software  \n";
	print $fh "  muscle = /home/software/muscle3.8.31/muscle3.8.31\n";
	print $fh "\n";
	print $fh "# genome files TO EXTRACT SEQUENCES (can be symbolic links created with ln -s)\n";
	print $fh "  #Note that you can simply create symbolic links of the genome files, with ln -s\n";
	print $fh "  genone_a   =  /data/DelGet/MyRun/genomes_toalign/genome1.fa\n";
	print $fh "  gentwo_a   =  /data/DelGet/MyRun/genomes_toalign/genome2.fa\n";
	print $fh "  genthree_a =  /data/DelGet/MyRun/genomes_toalign/genome3.fa\n";
	print $fh "  \n";
	print $fh "# [FYI]: alignment software will be decided based on anch_dist\n";
	print $fh "  #ALNSOFT = muscle if anch_dist < 40000\n";
	print $fh "  #ALNSOFT = kalign if anch_dist >= 40000\n";
	print $fh "  \n";
	print $fh "##########################################################\n";
	print $fh "# FOR MASKING\n";
	print $fh "##########################################################\n";
	print $fh "# repeat masker software\n";
	print $fh "  RMSOFT = /home/software/RepeatMasker/RepeatMasker\n";
	print $fh "\n";
	print $fh "# masking librarie(s)\n";
	print $fh "# If blank, no masking will be done. \n";
	print $fh "# Several libraries can be defined, just separate by a , in the mask and the masking_folders list\n";
	print $fh "# As many masking_folders elements as mask files are required\n";
	print $fh "# -lib is needed in running repeat masker => use Repeat Masker script in their utilities to convert their embl to fasta if you would like to mask with the whole library\n";
	print $fh "  mask = (/data/DelGet/TEdb_to_mask/RepeatMaskerLib.20140131.fa,/data/DelGet/TEdb_to_mask/TE.subset.fa)\n";
	print $fh "  masking_folders = (masked.RM20140131,masked.subset)\n";
	print $fh "\n";
	close $fh;
	
	print STDERR "\n   Configuration file example printed in $cfile - once editited, run: perl DelGet.pl DelGet--0--CONFIG_manipname.pl\n\n";
	
	return 1;
}

#----------------------------------------------------------------------------
# print config file if asked
# load_config($conf);
#----------------------------------------------------------------------------
sub load_config {
	my $cfile = shift;	
	my %cf = ();
	open(my $fh, "<", $cfile) or die "\t    ERROR: SUB print_config: could not open to write $cfile!\n";
	CONFIG: while(<$fh>) {
		chomp (my $line = $_);
		$line =~ s/\s//g;
		next CONFIG if ((substr($line,0,1) eq "#") || ($line !~ /\w/) || ($line !~ /=/));
		my @line = split("=",$line);
		$line[1] = "" unless ($line[1]);
		$line[1] =~ s/#.*$//;
		$cf{$line[0]}=$line[1];
	}
	#Now the ones that are calculated:
	($cf{'anch_dist'} >= 40000)?($cf{'ALNSOFT'} = $cf{'kalign'}):($cf{'ALNSOFT'} = $cf{'muscle'});
	$cf{'a_max_len'} = $cf{'anch_len'} + (50 * $cf{'anch_len'} /100);
	$cf{'minlen'} = $cf{'anch_dist'} + (2 * $cf{'anch_len'});	
	$cf{'multip'} = 7 if ($cf{'anch_dist'} <= 10000); #10kb intially => 70kb max [muscle]\n";
	$cf{'multip'} = 4 if (($cf{'anch_dist'} <= 50000) && ($cf{'anch_dist'} > 10000)); #25kb => 100kb max [muscle], 50kb => 200kb max [kalign software]
	
	close $fh;	
	return(\%cf);
}

#----------------------------------------------------------------------------
# get assembly gaps coordinates, from standard UCSC format, see below
# $gapgen_ref = DelGet::get_gaps($gapgen_ref,$file,$gapfiles[$f]); #access to gaps of $file = $gapgen_ref->{$file} since it's refs
#----------------------------------------------------------------------------
#0		1							2			3			4	5	6
#bin	chrom						chromStart	chromEnd	ix	n	size	type	bridge
#.		gi|432120764|gb|KB097841.1|	268			451			.	.	184		.		.
sub get_gaps {
	my ($gapgen_ref,$ID,$gapfile) = @_; #genome numer (ordered in a list from config file)
	open(my $fh, "<", $gapfile) or die "\t    ERROR: SUB get_gaps: could not open $gapfile!\n";
	GAPGEN: while(<$fh>) {
		chomp (my $line = $_);
		next GAPGEN if (($line =~ /#bin/) || ($line !~ /\w/));
		my @gapline = split(/\t/,$line);
		my $chr = $gapline[1];
		(exists $gapgen_ref->{$ID}{$chr})?($gapgen_ref->{$ID}{$chr} = "$gapgen_ref->{$ID}{$chr},$gapline[2],$gapline[3]"):($gapgen_ref->{$ID}{$chr} = "$gapline[2],$gapline[3]");
	}
	close $fh;
	return $gapgen_ref;
}

#----------------------------------------------------------------------------
# get previous regions
# => $r = get_previous_OKreg($OKregions,$log,$c,$prev_reg);
#----------------------------------------------------------------------------
# Memorize regions coordinates for the ones in genome1, ie where anchors are picked
# Structure of _OKregions.tab output file starts like that:
# 0		1		2		3		4		5		6		7		8		9
# ID	type	Gname	Gstart	Gend	type	Gname	Gstart	Gend	DIST
sub get_previous_OKreg {
	($path,$r,$OKregions,$log,$c,$prev_reg) = @_;
	open(my $log_fh, ">>", $log) or print "ERROR: could not open to write $log $!\n";
	print $log_fh "\t--- Previous outputs are set to be loaded to avoid overlaps (RUN $c):\n";
	if ($c == 0) { #meaning no previous Deletions.X folders
		if ($OKregions eq "no"){ #no previous file
			print $log_fh "\t    => no file will be loaded to avoid overlaps (because this is first run and no previous file path defined in config file)\n";
		} else { #previous file defined => will be loaded
			print $log_fh "\t    => $OKregions will be loaded to avoid overlaps (previous file path defined in config file)\n";
		}
	} elsif ($c == 1) { #Second run OR there was a Deletions.0 folder already; in both cases, get the previous OK regions
		if ($OKregions eq "no"){ #no previous file => OKregions in Deletions.0 is $OKregions
			$OKregions = "$path/Deletions.0/_OKregions.tab";
		} else { #need to cat the defined $OKregions in config file with OKregions in Deletions.0 + replace value of variable to load the correct file
			system "cat $path/Deletions.0/_OKregions.tab $OKregions > $path/Deletions.0/_OKregions.all.tab";
			$OKregions = "$path/Deletions.0/_OKregions.all.tab";
		}
	} elsif ($c == 2){
		my $check = $OKregions;
		$OKregions = "$path/Deletions.1/_OKregions.all.tab";
		system "cat $path/Deletions.0/_OKregions.tab $path/Deletions.1/_OKregions.tab > $OKregions" if ($check eq "no");
		system "cat $path/Deletions.0/_OKregions.all.tab $path/Deletions.1/_OKregions.tab > $OKregions" if ($check ne "no");
	} elsif ($c != 0) { #more than 2 => cat prev .all.tab and new one
		my $c_prev1 = $c - 1;
		my $c_prev2 = $c - 2;
		my $prev1path = "$path/Deletions.$c_prev1";
		my $prev2path = "$path/Deletions.$c_prev2";
		$OKregions = "$prev1path/_OKregions.all.tab";
		system "cat $prev2path/_OKregions.all.tab $prev1path/_OKregions.tab > $OKregions";
	}
	print $log_fh "\t    => $OKregions will be loaded to avoid overlaps\n" unless ($c == 0);
	
	# load OK region file unless not relevant, ie no previous runs
	if (($c != 0) && ($OKregions ne "no")) {#besides this situation, there will be a file to load
		open(my $prev_reg_fh, ">", $prev_reg) or die print $log_fh "\t    ERROR - can not create file $prev_reg $!";
		print $prev_reg_fh "#chr\tstart\tend\n\n";
		open(my $OKregions_fh, "<", $OKregions) or die print $log_fh "\t    ERROR - can not open previous _OKregions output file $OKregions $!";
		my $reg_r;
		PREVREG: while(<$OKregions_fh>){
			chomp (my $line = $_);
			next PREVREG if (($line =~ /^#/) || ($line !~ /\w/));
			my @line = split(/\t/,$line);
			print $prev_reg_fh "$line[2]\t$line[8]\n";
			
			# get the round to start with
			my $curr_reg_round = $line[0];
			$curr_reg_round =~ s/^reg(.*)-.*$/$1/; #extract this line round number
			$reg_r = $curr_reg_round unless ($reg_r); #first round, has to be remembered
			$reg_r = $curr_reg_round if ($reg_r < $curr_reg_round); #memorize region ID, if larger (that way still OK even if not numerical order in the file, i.e. excel sorts by alphabetical)
		}
		$r = $reg_r + 1; #reinitialize $r to this last round number +1 => still unique IDs for regions
		print $log_fh "\t\t=> first round of this run = $r (ie last round of the previous run was $reg_r) => new regions will start at $r\n";
	}
	close $log_fh;
	return $r;
}

#----------------------------------------------------------------------------
# load previous regions
# => my $alreadyrand_full_ref = DelGet::load_previous_OKreg($prev_reg);
#----------------------------------------------------------------------------
sub load_previous_OKreg {
	my $prev_reg = shift;
	my %alreadyrand_full = ();
	open(my $prev_reg_fh, "<", $prev_reg) or die "\t    ERROR - can not open file $prev_reg $!";
	PREVREG_IN: while (<$prev_reg_fh>) {
		chomp (my $line = $_);
		next PREVREG_IN if (($line =~ /^#/) || ($line !~ /\w/));
		my @line = split(/\t/,$line);
		(exists $alreadyrand_full{$line[0]})?($alreadyrand_full{$line[0]} .= ",$line[1]"):($alreadyrand_full{$line[0]} = $line[1])
	}
	close $prev_reg_fh;
	return \%alreadyrand_full;
}	

#----------------------------------------------------------------------------
# Check for overlap between regions
# ($ifnextrand) = DelGet::check_overlap_reg($randomized_ends{$rdmID},$anch_len,$anch_dist,$three_end,$five_start) if (exists $randomized_ends{$rdmID});
#----------------------------------------------------------------------------
sub check_overlap_reg {
	my ($already,$anch_len,$anch_dist,$three_end,$five_start) = @_;
	my $ifnext = "no";
	my @already = split(",",$already);
	REGCHECK: foreach my $prev_end (@already) {
		my $prev_start = $prev_end - (2*$anch_len) - $anch_dist + 3;
		if ((($three_end < $prev_end) && ($three_end > $prev_start)) || (($five_start < $prev_end) && ($five_start > $prev_start))) {
			$ifnext = "yes";
			last REGCHECK;
		}
	}
	return ($ifnext);
}

#----------------------------------------------------------------------------
# Check for overlap between randomized region and assembly gap
# ($ifnextrand) = DelGet::check_overlap_gap($gapgen_ref->{0}->{$rdmID},$three_end,$five_start) if (exists $gapgen_ref->{0}->{$rdmID});
# my ($ifnext2) = DelGet::check_overlap_gap($gapgen_ref->{1}->{$t_name_gen2},$end_gen2,$start_gen2);
# my ($ifnext3) = DelGet::check_overlap_gap($gapgen_ref->{2}->{$t_name_gen3},$end_gen3,$start_gen3);
#----------------------------------------------------------------------------
sub check_overlap_gap {
	my ($storedgaps,$three_end,$five_start) = @_;
	my $ifnext = "no";
	my @gaps = split(",",$storedgaps);
	GAPCHECK: for (my $g=0;$g<=$#gaps;$g=$g+2) {
		my $gap_start = $gaps[$g];
		my $gap_end = $gaps[$g+1];
		if (($gap_end > $five_start) && ($gap_start < $three_end)) {
			$ifnext = "yes"; #there is overlap, no need to check the rest => will skip this region
			last GAPCHECK;
		}	
	}
	return ($ifnext);
}

#----------------------------------------------------------------------------
# Extract sequences + get infos of anchor coords for final output
# => my $gen1Infos_ref = DelGet::extract_sequences($tempposi,$anchors,$log,$genone);
#----------------------------------------------------------------------------
sub extract_sequences {
	my ($tempposi,$anchors,$log,$genome) = @_;
	open(my $log_fh, ">>", $log) or print "ERROR: could not open to write $log $!\n";	
	open(my $tempposi_fh,"<", $tempposi) or die print $log_fh "\t    ERROR - can not open file $tempposi $!";
	my $anchorsseqio = Bio::SeqIO->newFh(-format => 'Fasta', -file=>">$anchors") or die print $log_fh "Failed to create SeqIO FH object from $anchors $!\n";
	my $db = Bio::DB::Fasta->new($genome) or die print $log_fh "\t\tERROR: could not create Bio::DB::Fasta object from $genome $!\n";
	print $log_fh "\t--- extracting sequences...\n";
	my %gen1Infos = ();
	EXT: while(<$tempposi_fh>){
		next EXT unless ($_ =~ /\w/);	#avoid blank lines
		chomp (my $line = $_);	
		my ($anchor_name,$Gname,$Gstart,$Gend)= split(/\t/,$line);
		my $subSeq = $db->seq($Gname,$Gstart,$Gend) or print $log_fh "\t\t$Gname was not found in $genone\n";	# extract target sequence, unless problem
		my $seqobj = Bio::Seq->new( -display_id => $anchor_name, -seq => $subSeq); # create object with target
		print $anchorsseqio $seqobj;	# print it out (in fasta format)
		
		#keep coords in memory for final output
		$gen1Infos{$anchor_name} = "$Gname\t$Gstart\t$Gend";
	}
	close $tempposi_fh;
	close $log_fh;
	return \%gen1Infos;
}	

#----------------------------------------------------------------------------
# parse blat outputs
# => my @HSfiles_temp = DelGet::parse_blat($log,$blatfiles[$i],$genIDs[$i],$pathtemp);
#----------------------------------------------------------------------------
sub parse_blat {
	my ($log,$file,$genID,$pathtemp) = @_;
	open(my $log_fh, ">>", $log) or print "ERROR: could not open to write $log $!\n";	
	print $log_fh "\t    parsing $file in progress...\n";
	
	#going through input + writing output
	my %highestscore = ();
	my %howhighest = ();
	my %blatline = ();
	my $lc = 0;
	open($file_fh,"<",$file) or die print $log_fh "\t    ERROR - can not open blat output file $file $!";
	while(<$file_fh>){
		chomp (my $line = $_);
		unless ($lc < 5) {
			my @line = split(/\s+/,$line);
			my ($matches,$mismatches,$rep_matches,$q_name,$q_length) = ($line[0],$line[1],$line[2],$line[9],$line[10]);
			my $score = sprintf "%.2f", ( 100 * ( $matches + $mismatches + $rep_matches ) / $q_length );
			my ($region,$type) = split("#",$q_name);
			my $ID = $region."#".$genID."#".$type;
			if (exists $highestscore{$ID}) {
				if ($score > $highestscore{$ID}) {
					$howhighest{$ID} = $highestscore{$ID} / $score;
					$highestscore{$ID} = $score;
					$blatline{$ID} = $line;
				}
			} else { #first line basically
				$highestscore{$ID} = $score;
				$howhighest{$ID} = 1;
				$blatline{$ID} = $line;
			}
		}
		$lc ++;
	}
	close $file_fh;
	#write an output to keep track
	my @HSfiles_temp = ();
	foreach my $scID (sort keys %highestscore) {
		my ($region,$genID,$type) = split("#",$scID);
		my $OnlyHighScores = "$pathtemp/data_highest-scores/$region.HighestScores.tab";
		push (@HSfiles_temp,$OnlyHighScores);
		open($OnlyHighScore_fh,">>",$OnlyHighScores) or die print $log_fh "\t    ERROR - can not create output file $OnlyHighScores $!";
		print $OnlyHighScore_fh "$blatline{$scID}\t$region\t$type\t$genID\tscore=\t$highestscore{$scID}\tratio_with_2nd=\t$howhighest{$scID}\n";
	}
	print $log_fh "\t\t-> written in $pathtemp/data_highest-scores/regXX.HighestScores.tab...\n";
	close $log_fh;
	return (@HSfiles_temp);
}

#----------------------------------------------------------------------------
# get distance between the anchors
# => my ($undef,$dist_gen2) = DelGet::check_anchor_dist($strand_gen2,$gentwo_ref);
#----------------------------------------------------------------------------
sub check_anchor_dist { 
	($strand,%gen) = @_;
	my $undef = 0;
	my ($start,$end,$dist);
	if ($strand eq "+") {
		$start = $gen{5}->[2];
		$end = $gen{3}->[3];
		$dist = $gen{3}->[2] - $gen{5}->[3] + 1;
	} elsif ($strand eq "-") {
		$start = $gen{3}->[2];
		$end = $gen{5}->[3];
		$dist = $gen{5}->[2] - $gen{3}->[3] + 1;
	} else {
		$undef = 1; 
	}
	return ($undef,$start,$end,$dist);
}

	
	
#ensure last returned value is true (to load as require in scripts)
1;