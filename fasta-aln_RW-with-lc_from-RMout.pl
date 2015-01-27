#!/usr/bin/perl -w

#######################################################
# Author :  Aurelie K
# date   :  May 2013    
# email  :  4urelie.k@gmail.com
# Pupose :  Use an aligned file and a Repeat Masker output file (of the same sequences that were aligned) to rewrite sequences masked in lowercases
#####################################################
use warnings;
use Bio::Perl;
use Bio::SeqIO;

my $usage = "
	perl <scriptname.pl> <path>
	
	Typically, path is a directory were are located all .masked and .out files to be processed
	Note that only masked files and matching .out files will be [if file had no TEs in its sequences, won't be considered]\n";
	
my $path = shift @ARGV or die "$usage";

print "\nScript started\n";

#########################################################################################################
# Get list of regions to deal with, ie the ones that have a .masked output
##########################################################################################################
my @allfiles = `ls $path`;
my @maskedfiles = grep(/.*.masked$/,@allfiles);
my @maskedreg = ();
for (my $i = 0; $i <= $#maskedfiles; $i++) {
	$maskedfiles[$i] =~ s/\n//;
	my $reg = $maskedfiles[$i];
	$reg =~ s/^(reg[0-9]+-[0-9]+).fa.*/$1/;
	push(@maskedreg,$reg);
}


#########################################################################################################
# Now loop on these regions
##########################################################################################################
foreach my $region (@maskedreg) {
	# open files
	my $aln = "$path/$region.fa.align.fa";
	my $aln_in = Bio::SeqIO->new(-format => 'Fasta', -file=>$aln) or die "Failed to create Bio::SeqIO object from $aln $!\n";
	print " --- Dealing with region = $region ($aln)\n";
	#####################################################
	# get the gaps for each sequence + their length => being able to modify coordinates
	#####################################################
	# 123456789
	# --GTTTATCA => gap from 1 to 3, len = 2
	# TAGTT---CA => gap from 6 to 9, len = 3
	# TAGTTTTTCA
	#
	print "     Get gaps...\n";
	my %gaps = ();
	while( my $seq = $aln_in->next_seq() ) {
		my $header = $seq->display_id;
		my $alnlen = $seq->length();
		my @seq = split(//,$seq->seq);
		my $start;
		for (my $n=1;$n<$alnlen;$n++) {
			#gap opening: get start
			$start = $n if (($n == 1) && ($seq[$n-1] eq "-")); #just case of first nt of aln
			$start = $n+1 if (($seq[$n] eq "-") && ($seq[$n-1] ne "-"));
			if (($seq[$n-1] eq "-") && ($seq[$n] ne "-")) { 
			#gap ending: get end + print coords
				my $end = $n+1;
				my $len = $end - $start;
				$gaps{$header}.="$start\t$len\t$end\t";
			}
		}
	}
	
	#####################################################
	# deal with RMout, add pairs of coords that need to be rewritten
	#####################################################
	# 1 2 3 4 5 6 7 8 9 10
	# A C g t t t A T C A => masked from 3 to 6
	#     x x x x	
	my %hash_of_posi = ();
	my ($name,$start,$stop);
	my ($prevname,$prevstart,$prevstop) = ("firstline","firstline","firstline");
	my $count = 1;
	my $rmout = "$path/$region.fa.out";
	chomp (my $nb_lines = `grep -c "" $rmout`);
	print "     Get RMout coordinates (in $rmout - $nb_lines lines...\n";
	open(RMOUT,"<$rmout") or die "can not open file $rmout $!\n";
	while(<RMOUT>){
		chomp (my $line = $_);
		unless (($line =~ /.*query|sequence.*/) || ($line !~ /\w/) ){ #avoid blank lines + columns labels
			$line =~ s/^\s+//;					#To remove space(s) in the begining of each line (otherwise, problem with split)
			my @col = split(/\s+|\t/, $line);				#each element of the line separated by spaces are put in a list	
			($name,$start,$stop) = ($col[4],$col[5],$col[6]);
			if ($prevname eq "firstline") { #first line of this name
				$hash_of_posi{$name}.="$start\t";
			} elsif ($name eq $prevname) {
				unless ($start <= $prevstop) { #ie if overlap => nothing is written
					$hash_of_posi{$name}.="$prevstop\t$start\t";
				}	
			} else { #ie change of sequence => I need to add the stop to prevname hash
				$hash_of_posi{$prevname}.="$prevstop\t";
				$hash_of_posi{$name}.="$start\t";
			}
			if ($count == $nb_lines) { #ie last line of file
				$hash_of_posi{$name}.="$stop\t";
			}
			($prevname,$prevstart,$prevstop) = ($name,$start,$stop); #set values for next comparison
		}
		$count++;
	}
	close RMOUT;
	
	
	#####################################################
	# Modify RMout coords to add gaps
	#####################################################
	print "     Modify RMout coordinates...\n";
	my %posi_with_gaps = ();
	foreach my $seqname (sort keys %hash_of_posi) {
		my @masked_coords = split (/\t/,$hash_of_posi{$seqname});
		# in rare cases of no gaps in the sequence
		if ($gaps{$seqname}) {
			my @gaps_coords = split (/\t/,$gaps{$seqname});
			for (my $i=0; $i<=$#masked_coords; $i+=2) {
				my $new_start = $masked_coords[$i];
				my $new_end = $masked_coords[$i+1];
				for (my $j=0; $j<=$#gaps_coords; $j+=3) {
					$new_start += $gaps_coords[$j+1] if ($new_start > $gaps_coords[$j]);
					$new_end += $gaps_coords[$j+1] if ($new_end > $gaps_coords[$j]);
				}
				#re write positions; gap length has been added if needed
				$posi_with_gaps{$seqname}.="$new_start\t$new_end\t";
			}	
		} else {
			$posi_with_gaps{$seqname}=$hash_of_posi{$seqname};
		}
	}
	
	#KEEPING TRACK
	my $log = "$aln.masked.log";
	open(LOG, ">$log") or die "Failed to create file $log $!\n";
	print LOG "\nNew coordinates for masked regions = \n";
	foreach my $seqname (sort keys %posi_with_gaps) {
		print LOG "$seqname\t$posi_with_gaps{$seqname}\n";
	}
	print LOG "\n\nExtraction log (modified sequences + errors)\n";
	
	#####################################################
	# deal with aligned file, and rewrite.
	#####################################################
	# open files
	$aln_in = Bio::SeqIO->new(-format => 'Fasta', -file=>$aln) or die "Failed to create Bio::SeqIO object from $aln $!\n";
	open(OUT, ">$aln.masked") or die "Failed to create file $aln.masked $!\n";
	#deal with sequences in the alignement one by one	
	print "     Rewriting sequences...\n";
	while( my $seq = $aln_in->next_seq() ) {			
		my $header = $seq->display_id;
		my $len = $seq->length();
		print OUT ">$header\n";
		print LOG "$header\t";
		if (exists $posi_with_gaps{$header}) {
			my @coords = split (/\t/,$posi_with_gaps{$header});
	#		print "\n\nlength of $header is $len\n";
			unless ($coords[0] == 1) { #ie unless no first part to write in uppercase
				my $uc_first_stop = $coords[0]-1;
				my $firstpart = uc($seq->subseq(1,$uc_first_stop));
	#			print "firstpart; de 1 à $uc_first_stop, $firstpart \n";
				print OUT "$firstpart";
			}
			my $j = $#coords-2; #get number of elements -2 to loop on that
			for (my $i = 0; $i<=$j; $i+=2 ) { # loop on number of characters in the sequence
				my $lc_start = $coords[$i];
				my $lc_stop = $coords[$i+1];
				my $uc_start = $coords[$i+1]+1; #otherwise rewrite letter
				my $uc_stop = $coords[$i+2]-1; #otherwise rewrite letter
	#			print "THIS LOOP TURN = $lc_start, $lc_stop, $uc_start, $uc_stop\n";
				if ($lc_stop == $uc_stop) { #ie the two lower cases are adjacent => no upercase
					unless ($lc_start > $lc_stop) {
						my $lc_seq = lc($seq->subseq($lc_start,$lc_stop));
	#					print "$lc_stop == $uc_stop => $lc_start to $lc_stop = $lc_seq\n";
						print OUT "$lc_seq";
					} else {
						print LOG "$lc_start > $lc_stop\n";
					}
				}
				if ($lc_stop == $coords[$i+2]) { #ie prob one letter would be written too times => need to be modified
					$lc_stop = $lc_stop-1; #correct that
					unless ($lc_start > $lc_stop) {
						my $lc_seq = lc($seq->subseq($lc_start,$lc_stop));
	#					print "$lc_stop = $coords[$i+2] => $lc_start to $lc_stop = $lc_seq\n";
						print OUT "$lc_seq";
					} else {
						print LOG "$lc_start > $lc_stop\n";
					}	
				} elsif ($lc_stop != $uc_stop) {	
					unless ($lc_start > $lc_stop) {
						my $lc_seq = lc($seq->subseq($lc_start,$lc_stop));
	#					print "else lc - $lc_start to $lc_stop = $lc_seq\n";
						print OUT "$lc_seq";
					} else {
						print LOG "$lc_start > $lc_stop\n";
					}	
					unless ($uc_start > $uc_stop) {
						my $uc_seq = uc($seq->subseq($uc_start,$uc_stop));
	#					print "else uc - $uc_start to $uc_stop = $uc_seq\n";
						print OUT "$uc_seq";
					} else {
						print LOG "$uc_start > $uc_stop\n";
					}	
				}	
			}
			my $lc_last_start = $coords[-2];
			my $lc_last_stop = $coords[-1];
			my $lc_last = lc($seq->subseq($lc_last_start,$lc_last_stop));
	#		print "$lc_last_start to $lc_last_stop = $lc_last\n";
			print OUT "$lc_last";
			unless ($lc_last_stop == $len) { #ie no last part to write
				my $uc_last_start = $coords[-1]+1;
				my $uc_last = uc($seq->subseq($uc_last_start,$len));
	#			print "lastpart; de $uc_last_start à $len, $uc_last \n";
				print OUT "$uc_last";
			}	
			print OUT "\n";
			print LOG "\n";
		} else {
			my $unmasked = $seq->seq();
			print OUT "$unmasked\n";
		}
	}
	close OUT;

	# Rewrite genome in fasta format
	my $output = "$aln.masked.fa";
	my $maskedaln = Bio::SeqIO->new(-file => "<$aln.masked", -format => "fasta") or die "Failed to create Bio::SeqIO object from $aln.masked $!\n";
	my $maskedalnfasta = Bio::SeqIO->new(-file => ">$output", -format => "fasta") or die "Failed to create Bio::SeqIO outputfile $output $!\n";
	while( my $seq = $maskedaln->next_seq() ) {
		$maskedalnfasta->write_seq($seq);		
	}
	unlink "$aln.masked";
}
print "\nFinished, see output files masked.fa in $path.
In .log files are listed modified sequences.\n\n";
exit;






