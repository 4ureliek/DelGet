#!/usr/bin/perl -w

##########################################################
# Author  :  Aurelie K
# version :  1.0 (see updates)
# email   :  4urelie.k@gmail.com
# PURPOSE :  After running pipeline, get all gap freq outputs concatenated 
##########################################################
# UPDATES
#	- v1.0 = 17 Oct 2013
#	- v1.1 = 05 Nov 2013
#		Error in getting file name for lcuc parsing concatenation => was keeping only last folder
#	- v1.2 = 03 Mar 2014
#		problem when no gaps in some categories when concat.
##########################################################
use strict;
use warnings;


#################################################################
# Usage and get argument
#################################################################
my $usage = "\nUSAGE:
	perl script.pl path
	
	Where path = folder where Deletions.X directories are, corresponding to a single 'experiment'
	Without the \"/\" at the end of the path

	Typically, if a path has to be defined:
		perl DelGet.pl /data/XXX/Deletion_pipeline
	OR, if current directory:
		perl DelGet.pl .\n\n";
	

my $path = shift @ARGV or die "$usage $!\n" ;

print "\n --- script started\n";
#################################################################
# Get folders where to look for
#################################################################
my @EAdir = `ls -d $path/Deletions.*/_ExtractAlign/` or die print " !! ERROR - can't list files in $path $!\n"; #Deletions.X/_ExtractAlign folders only, full path

#################################################################
# Gather the gapfreq(.lcuc).log files and store stuff to concat
#################################################################
# have coordinates for different values, to append to "type" => print in good order
my %coords = (
	'Total length of alignements' => '0',
	'- Total number of spe gaps 1-30nt' => '1',
	'- Total length of spe gaps 1-30nt' => '2',
	'- Total number of spe gaps >30nt' => '3',
	'- Total length of spe gaps >30nt' => '4',
);

# loop on Deletions.X folders
my %gapfreq = ();
my %lcuc = ();
my @RMout = ();
my $i = 0;
print "\n --- Processing gapfreq results...\n";
foreach my $EAdir (@EAdir) {
	chomp $EAdir;
	$EAdir =~ s/\/$//;
	print "     -> $EAdir\n";
	#If first folder
	if ($i == 0){
		#get the 2 species
		my @files =`ls $EAdir`;
		foreach my $file (@files) {
			chomp (my $spec = $file);
			$spec =~ s/^_([A-Z][a-z][a-z][a-z])\.gaps\.specific\.bed$/$1/;
			$gapfreq{$spec}=1 if ($file =~ /^_[A-Zaz]+\.gaps\.specific\.bed$/);
		}
		$i=1;
	}
	#get masked folders if they exist
	my @RMout_temp = `ls -d $EAdir/*/` or die print " ERROR - can't list files in $EAdir $!\n"; #get folders with masked stuff
	foreach my $masked (@RMout_temp) {
		$masked =~ s/\/$//;
		push(@RMout,$masked) if ("$masked/gapfreq.lcuc.log"); #push only if masked files there that were analyzed
	}
	
	#there should be a gapfreq.log in this folder
	open(GAPFREQ,"<$EAdir/gapfreq.log") or die print " !! ERROR: could not open file $EAdir/gapfreq.log $!\n";
	my $file;
	my $parse = 0;
	GF: while (<GAPFREQ>) {
		chomp (my $temp = $_);
	
		#skip anything that's not results
		$parse = 1 if ($temp =~ /.*---\sparsing\sgaps.*/);
		next GF unless ($parse == 1);
		
		#get file if relevant
		$file = $temp if ($temp =~ /_[A-Za-z]+\.gaps\.specific\.bed/);
		$file =~ s/.*\/_ExtractAlign\/(_[A-Za-z]+\.gaps\.specific\.bed)/$1/ if ($temp =~ /_[A-Za-z]+\.gaps\.specific\.bed/);

		#next should be values, split on tab
		if (($temp =~ /.*Total.*of\sspe\sgaps.*/) || ($temp =~ /.*Total\slength\sof\salignements.*/)){
			my ($type,$value) = split(/\t/,$temp);
			$type =~ s/^\s+//;
			$type = $coords{$type}." ".$type;
			(defined $gapfreq{$file}{$type})?($gapfreq{$file}{$type}+=$value):($gapfreq{$file}{$type}=$value);
		}
		#when that's done, it will be next species or end of file
	}
	close GAPFREQ;
}#end loop Deletions.X folder for gapfreq


#now if relevant also get the lcuc
print "\n --- Processing gapfreq.lcuc results if any...\n";
foreach my $masked (@RMout) {
	chomp $masked;
	my $file_lcuc;
	my $parse_lcuc = 0;
	print "     -> $masked\n";
	open(GAPFREQ_LCUC,"<$masked/gapfreq.lcuc.log") or die print " !! ERROR: could not open file $masked/gapfreq.lcuc.log $!\n";
	GF_LCUC: while (<GAPFREQ_LCUC>) {
		#skip anything that's not results
		chomp (my $temp_lcuc = $_);
		$parse_lcuc = 1 if $temp_lcuc =~ /--- parsing gaps/;
		next GF_LCUC unless ($parse_lcuc == 1);
		
		#get file if relevant
		if ($temp_lcuc =~ /_[A-Za-z]+\.gaps.*\masked.*bed/) {
			$file_lcuc = $temp_lcuc;
			$file_lcuc =~ s/.*\/_ExtractAlign\/.*masked.*\/(_[A-Za-z]+\.gaps\.in.*masked.*bed)/$1/;
			$file_lcuc =~ s/masked\.bed/masked\.specific\.bed/ unless $file_lcuc =~ /specific/;
		}
		#next should be values, split on tab
		if (($temp_lcuc =~ /Total.*of\sspe\sgaps/) || ($temp_lcuc =~ /Total\slength\sof\salignements/)){
			my ($type,$value) = split(/\t/,$temp_lcuc);
			$type =~ s/^\s+//;
			$type = $coords{$type}." ".$type;
			my $masked_dirname = $1 if ($masked =~ m/.*Deletions\..*\/_ExtractAlign\/(.*masked.*)/);
			(defined $lcuc{$masked_dirname}{$file_lcuc}{$type})?($lcuc{$masked_dirname}{$file_lcuc}{$type}+=$value):($lcuc{$masked_dirname}{$file_lcuc}{$type}=$value);
		}
	}
	close GAPFREQ_LCUC;
}#end loop foreach RMout (masked folders -> lcuc gaps)


print "\n --- done getting values\n\n";

#################################################################
# Now read hashes and print cat results
#################################################################
print " --- Reading hashes and printing outputs\n";
#gapfreq
my $out = "$path/_cat.gapfreq.out";
open(OUT,">$out") or die print " !! ERROR: could not create file $out $!\n";
foreach my $file (keys %gapfreq) {
	print OUT "$file\n";
	my $t;
	foreach $t (sort keys %{$gapfreq{$file}}) {
		print OUT "$t\t$gapfreq{$file}{$t}\n";
	}
}
close OUT;

#lcuc files
foreach my $masked (@RMout) {
	my $masked_dirname = $1 if ($masked =~ m/.*Deletions\..*\/_ExtractAlign\/(.*masked.*)/);
	my $out_lcuc = "$path/_cat.gapfreq.$masked_dirname.lcuc.out";
	open(OUT_LCUC,">$out_lcuc") or die print " !! ERROR: could not create file $out_lcuc $!\n";
	foreach my $file_lcuc (keys %{$lcuc{$masked_dirname}}) {
		print OUT_LCUC "$file_lcuc\n";
		my $t;
		foreach $t (sort keys %{$lcuc{$masked_dirname}{$file_lcuc}}) {
				print OUT_LCUC "$t\t$lcuc{$masked_dirname}{$file_lcuc}{$t}\n";
		}
	}
	close OUT_LCUC;
}
print " --- done [script done]\n\n\n";
exit;







