#!/usr/bin/perl -w

#######################################################
# SUBROUTINES
# FOR SUBSET OF SCRIPTS = DelGet.pl
# => defined as Deletions package
######################################################
package Deletions;

######################################################
# QUITE GENERAL FILE ARCHITECTURE STUFF
######################################################
#----------------------------------------------------------------------------
# get a filename from a full path
# my $genone_name = Deletions::filename($genone);
#----------------------------------------------------------------------------
sub filename {
	my($name) = shift;
	$name =~ s/.*\/(.*)$/$1/;
	return $name;
}
#----------------------------------------------------------------------------


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
#----------------------------------------------------------------------------


######################################################
# QUITE GENERAL FASTA STUFF
######################################################
#----------------------------------------------------------------------------
# calculate length of the genome, for sequences > $minlen - includes check steps to avoid repeating this if length already calculated
# => my ($totlength) = Deletions::get_tot_len_filtered($genome,$log,$minlen);
# 	 note that log file need to be closed in the main before calling the subroutine, and reopened after
#----------------------------------------------------------------------------
sub get_tot_len_filtered {
	my $genome = shift;
	my $log = shift;
	my $minlen = shift;
	
	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	

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
	} else {
		print $log_fh "\t\tTotal length of genome for scaffolds > $minlen has been previously calculated ($lengthfile exists)\n";
		open (my $len_fh2, "<", $lengthfile) or die print $log_fh "ERROR: could not open $lengthfile $!\n";
		while (<$len_fh2>) {
			$GenLen = $_;
		}	
	}
	print $log_fh "\t\t => total len = $GenLen nt\n";
	return($GenLen);
}
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Get all lengths of the genome - includes check steps to avoid repeating this if length already calculated
# => my ($geninfos) = Deletions::get_all_len_filtered($genome,$log,$minlen,$GenLen);
# 	 note that log file need to be closed in the main before calling the subroutine, and reopened after
#----------------------------------------------------------------------------
sub get_all_len_filtered {
	my ($genome,$log,$minlen,$GenLen,$lengthfile,$randtot) = @_;
	
	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	
	
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
				print $len_fh "$ID\t$len\t$ratio\n";
				my @infos = ($len,$ratio);
				$geninfos{$ID} = \@infos;
				push(@dbIDs_filtered,$ID);
			}
		}
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
	}
	return(\%geninfos,@dbIDs_filtered);
}
#----------------------------------------------------------------------------



######################################################
# QUITE SPECIFIC TO THIS PIPELINE
######################################################
#----------------------------------------------------------------------------
# get assembly gaps coordinates, from standard UCSC format, see below
# $gapgen_ref = Deletions::get_gaps($file,$gapfiles[$f]); #access to gaps of $file = $gapgen_ref->{$file} since it's refs
#----------------------------------------------------------------------------
#0		1							2			3			4	5	6
#bin	chrom						chromStart	chromEnd	ix	n	size	type	bridge
#.		gi|432120764|gb|KB097841.1|	268			451			.	.	184		.		.
sub get_gaps {
	my $ID = shift; #genome numer (ordered in a list from config file)
	my $input = shift; #file name
	my %gaps = ();
	open(my $input_fh, "<", $input) or die "\t    ERROR: SUB get_gaps: could not open $input!\n";
	GAPGEN: while(<$input_fh>) {
		chomp (my $line = $_);
		next GAPGEN if (($line =~ /#bin/) || ($line !~ /\w/));
		my @gapline = split(/\t/,$line);
		my $chr = $gapline[1];
		(exists $gaps{$ID}{$chr})?($gaps{$ID}{$chr} = "$gaps{$ID}{$chr},$gapline[2],$gapline[3]"):($gaps{$ID}{$chr} = "$gapline[2],$gapline[3]");
	}	
	return \%gaps;
}
#----------------------------------------------------------------------------


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
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";
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
	unless (($OKregions eq "no") && ($c == 0)) {#besides this situation, there will be a file to load
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
	return $r;
}
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# load previous regions
# => my $alreadyrand_full_ref = Deletions::load_previous_OKreg($prev_reg);
#----------------------------------------------------------------------------
sub load_previous_OKreg {
	my $prev_reg = shift;
	my %alreadyrand_full = ();
	if (-e $prev_reg) {
		open(my $prev_reg_fh, "<", $prev_reg) or die "\t    ERROR - can not open file $prev_reg $!";
		PREVREG_IN: while (<$prev_reg_fh>) {
			chomp (my $line = $_);
			next PREVREG_IN if (($line =~ /^#/) || ($line !~ /\w/));
			my @line = split(/\t/,$line);
			(exists $alreadyrand_full{$line[0]})?($alreadyrand_full{$line[0]} .= ",$line[1]"):($alreadyrand_full{$line[0]} = $line[1])
		}
	}	
	return \%alreadyrand_full;
}	
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Check for overlap between regions
# => my ($ifnextrand,$failed_anchor_ref) = Deletions::check_overlap_reg($randomized_ends{$rdmID},$anch_len,$anch_dist,$three_end,$five_start)
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
# => my ($ifnextrand,$failed_anchor_ref) = Deletions::check_overlap_gap($gapgen_ref->{1}->{$rdmID},$three_end,$five_start) if (exists $gapgen_ref->{1}->{$rdmID});
# => my ($ifnext2,$failed2) = Deletions::check_overlap_gap($gapgen_ref->{2}->{$t_name_gen2},$end_gen2,$start_gen2) if ($gapgen_ref->{2}->{$t_name_gen2});
#----------------------------------------------------------------------------
sub check_overlap_gap {
	my ($storedgaps,$three_end,$five_start) = @_;
	my $ifnext = "no";
	my @gaps = split(",",$storedgaps);
	GAPCHECK: for (my $g=0;$g<=$#gaps;$g=$g+2) {
		my $gap_start = $gaps[$g];
		my $gap_end = $gaps[$g+1];
		if (($gap_end > $five_start) && ($gap_start < $three_end)) {
			$ifnext = "yes";
			last GAPCHECK;
		}	
	}
	return ($ifnext);
}
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Extract sequences + get infos of anchor coords for final output
# => my $gen1Infos_ref = Deletions::extract_sequences($tempposi,$anchors,$log,$genone);
#----------------------------------------------------------------------------
sub extract_sequences {
	my ($tempposi,$anchors,$log,$genome) = @_;
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	
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
	return \%gen1Infos;
}	
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# parse blat outputs
# => my @HSfiles_temp = Deletions::parse_blat($log,$blatfiles[$i],$genIDs[$i],$pathtemp);
#----------------------------------------------------------------------------
sub parse_blat {
	my ($log,$file,$genID,$pathtemp) = @_;
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	
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
	return (@HSfiles_temp);
}
#----------------------------------------------------------------------------



#----------------------------------------------------------------------------
# get distance between the anchors
# => my ($undef,$dist_gen2) = Deletions::check_anchor_dist($strand_gen2,$gentwo_ref);
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