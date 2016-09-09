#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below / see change log
# email   :  4urelie.k@gmail.com
# PURPOSE :  DelGet util - to help checking the longest events
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;

my $version = "1.0";
my $scriptname = "DelGet--4--check_longest_events.pl";
my $changelog = "
#	- v1.0 = 07-08 Sept 2016
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname -i <directory> [-l] [-r] [-d] [-v] [-h] [-l]
    			
    PURPOSE:
    Find longest events, copy corresponding files, to help checking them
    Note that it expects that the configuration file be named .conf
		
    MANDATORY ARGUMENTS:	
     -i,--in      => (STRING) DelGet directory name (where the \"Del\" directory is located)
                              If several directory to parse, you can separate file names with commas: \",\"
                              or if -d is set, -i can designate a directory that contains the directories to parse 
 
    OPTIONAL ARGUMENTS
     -d,--dir     => (BOOL)   if -in is a directory containing the directories to parse
     -l,--lib     => (STRING) set a library to mask with (required, unless you set -s)
                              Default = /usr/local/RepeatMasker/Libraries/RepeatMaskerLib.20150807.fa
     -s,--sp      => (STRING) set species for -species in repeat masker 
                              instead of using the whole library above (to go faster)
     -r,--rm      => (STRING) repeat masker RepeatMasker script path
                              Default = /usr/local/RepeatMasker/RepeatMasker
     -u,--uc      => (STRING) Path to the script fasta-aln_RW-with-lc_from-RMout.pl in DelGet/Utils
                              Default = /home/akapusta/DelGet/Utils
     -v,--v       => (BOOL)   verbose, makes the script talk to you
     -v,--v       => (BOOL)   print version if only option
     -h,--help    => (BOOL)   print this usage
     -l,--log     => (BOOL)   print the change log (updates)
\n";

################################################################################
### Get arguments/options
### check some of them, print details of the run if verbose chosen
################################################################################
my $rm = "/usr/local/RepeatMasker/RepeatMasker";
my $lib = "/usr/local/RepeatMasker/Libraries/RepeatMaskerLib.20150807.fa";
my $sp = "na";
my $uc = "/home/akapusta/DelGet/Utils";
my ($in,$d,$help,$chlog,$v);
GetOptions ('in=s'      => \$in,  
            'dir'       => \$d, 
            'lib=s'     => \$lib, 
            'sp=s'      => \$sp, 
            'rm=s'      => \$rm,
            'uc=s'      => \$uc,
            'log'       => \$chlog, 
            'help'      => \$help, 
            'v'         => \$v);

#check step for options
die "\n $scriptname version: $version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n Please provide input directory (-i, use -h to see the usage)\n\n" if (! $in);
die "\n directory -i $in does not exist?\n\n" if (($in !~ /,/) && (! -e $in));
$in = $1 if ($in =~ /^(.*)\/$/); #remove the / at the end if any
die "\n Repeatmasker path incorrect?\n\n" if (! -e $rm);
die "\n Repeatmasker libraries path incorrect?\n\n" if (! -e $lib);
die "\n path to fasta-aln_RW-with-lc_from-RMout.pl incorrect?\n\n" if (! -e $uc);
$uc = $uc."/fasta-aln_RW-with-lc_from-RMout.pl";
#number of longest events to get
my $n = 5;

#output directory name (it will be created automatically in the same directory of $in)
my $checks = $in.".check_longest";

################################################################################
### MAIN
################################################################################
print STDERR "--------------------------------------------------\n" if ($v);
print STDERR " --- Get list of input director(y/ies)\n" if ($v);
($d)?($d="y"):($d="n");
#get list of input files + delete previous outputs
my ($listin,$outname) = input_files($in,$d);

# now loop through list
################################################################################
$checks = get_path($in)."/$checks";
print STDERR " --- All the alignment files will be copied in $checks\n" if ($v);
(-e $checks)?(`rm -Rf $checks`):(`mkdir $checks`);
print STDERR " --- Parsing files... \n" if ($v);
my $data = ();
INPUT: foreach my $in (@{$listin}) {
	chomp $in;
	print STDERR "      => From $in \n" if ($v);
	unless (-e $in) {
		print STDERR "         WARN: (main): $in does not exist? skipping\n";
		next INPUT;
	}
    my ($id1,$id2) = get_species($in);
	
	my $fext = "gaps.specific.bed.big.tab";
	my $id1_gaps = $in.".".$id1.".".$fext;
	print STDERR "         Concatenating $in/Del/Deletions.*/_ExtractAlign/_$id1.$fext in $id1_gaps \n" if ($v);
	`cat $in/Del/Deletions.*/_ExtractAlign/_$id1.$fext > $id1_gaps`;
	my $id2_gaps = $in.".".$id2.".".$fext;		
	print STDERR "         Concatenating $in/Del/Deletions.*/_ExtractAlign/_$id1.$fext in $id2_gaps \n" if ($v);
	`cat $in/Del/Deletions.*/_ExtractAlign/_$id2.$fext > $id2_gaps`;
	
	print STDERR "         Getting the $n longest events and copying alignment files\n" if ($v);
	$data = get_longest_del($in,$id1_gaps,$id1,$data,$n,$checks,$v);
	$data = get_longest_del($in,$id2_gaps,$id2,$data,$n,$checks,$v);

	print STDERR "         Printing the data in $outname.tab\n" if ($v);
	print_data($data,$outname,$v);
}	

print STDERR " --- Masking fasta files associated to aln files\n" if ($v);
mask_aln($checks,$rm,$lib,$sp,$uc,$v);
	
print STDERR "\n --- Done\n" if ($v);
print STDERR "--------------------------------------------------\n\n" if ($v);
exit;



################################################################################
### SUBROUTINES
################################################################################
#----------------------------------------------------------------------------
# from a filename or a directory keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#----------------------------------------------------------------------------
# get list of input files + delete previous outputs
# my ($listin,$outname) = input_files($in,$d);
#----------------------------------------------------------------------------
sub input_files {
	my($in,$d) = @_;
	my $path;
	my @listin = ();
	if ($d eq "y") {
		$path = $in;
		my @namelist = `ls $path`;
		NAME: foreach my $name (@namelist) {
			chomp($name);
			next NAME if ($name =~ /[Gg]enome/);
			$name = "$path/$name";
			push(@listin,$name) if (-d $name);
		}
		$outname = $path;
	} elsif ($in =~ /,/) {		
		@listin = split(",",$in);
		$path = get_path($listin[0]);
		$outname = $listin[0];
	} else {
		$outname = $in;
		push(@listin,$in);		
	}
	return (\@listin,$outname);
}

#----------------------------------------------------------------------------
# get species IDs
# my ($id1,$id2) = get_species($in);
#----------------------------------------------------------------------------
sub get_species {
	my($in,$v) = @_;
	my $conf = `ls $in/*.conf`;
	chomp($conf);
	if ($conf) {
		my $id1 = `grep IDgen2 $conf`;
		$id1 =~ s/\s//g;
		$id1 =~ s/#.*//;
		$id1 =~ s/IDgen2=//;
		my $id2 = `grep IDgen3 $conf`;
		$id2 =~ s/\s//g;
		$id2 =~ s/#.*//;
		$id2 =~ s/IDgen3=//;
		print STDERR "         The 2 species to consider = $id1,$id2\n" if ($v);
		return ($id1,$id2);
	} else {
		#conf file named something else...?		
		#to do: loop through files locted in there
		print STDERR "         WARN: (sub get_species): $in does not contain a .conf file - skipping\n";		
		return ("na","na");		
	}
}

#----------------------------------------------------------------------------
# get the longest deletion events
# $data = get_longest_del($in,$id1_gaps,$id1,$data,$n,$checks,$v);
#----------------------------------------------------------------------------
sub get_longest_del {
	my($in,$gaps,$id,$data,$n,$checks,$v) = @_;
	my @del = ();
	open(my $fh, "<", $gaps) or confess "     \nERROR (sub get_longest_del): could not open to read $gaps $!\n";
	#./Del/Deletions.0/_ExtractAlign/reg10-198.fa.align.fa  	6737   	6896   	67     	.      	+      	929    	10798
	#./Del/Deletions.0/_ExtractAlign/reg10-198.fa.align.fa  	7194   	7665   	67     	.      	+      	929    	10798
	while (<$fh>) {
		chomp($_);		
		my @l = split(/\s+/,$_);
		#need to re get the real length, after subtraction
		my $len = $l[2]-$l[1];
		my @info = ($l[0],$len);
		push(@del,\@info); 
	}
	close $fh;
	@del = sort { ($b->[1] <=> $a->[1]) } @del;		
	for (my $i=0; $i < $n; $i++) {
		my $aln = $del[$i]->[0];
		$aln =~ s/^\.\///;
		$aln = $in."/".$aln;
		$data->{$id}{$i}{$del[$i]->[1]}=$aln;
		my $aln_name = $checks."/".$id.".".$del[$i]->[1].".".filename($aln);
		`cp $aln $aln_name`;
		my ($fa,$fa_name) = ($aln,$aln_name);
		$fa =~ s/\.align\.fa$//;
		$fa_name =~ s/\.align\.fa$//;
		`cp $fa $fa_name`;
	}
	return($data);
}

#----------------------------------------------------------------------------
# print the longest events
# print_data($data,$outname,$v);
#----------------------------------------------------------------------------
sub print_data {
	my($data,$outname,$v) = @_;
	open(my $fh, ">", "$outname.tab") or confess "     \nERROR (sub print_json): could not open to write $outname.tab $!\n";
	foreach my $id (keys %{$data}) {
		foreach my $i (sort keys %{$data->{$id}}) {
			foreach my $len (sort keys %{$data->{$id}{$i}}) {
				print $fh "$id\t$i\t$len\t$data->{$id}{$i}{$len}\n";	
			}
		}
	}
	close $fh;
	return 1;
}

#----------------------------------------------------------------------------
# mask the alignments => bash script
# mask_aln($checks,$rm,$lib,$sp,$uc,$v);
#----------------------------------------------------------------------------
sub mask_aln {
	my ($checks,$rm,$lib,$sp,$uc,$v) = @_;
	my $cmd = "";
	($sp eq "na")?($cmd = "-lib $lib"):($cmd = "-sp $sp");
	my @files = `ls $checks`;	
	foreach my $f (@files) {
		chomp($f);
		$f = $checks."/".$f; #put back the path		
		
system( <<END_OF_BASH_SCRIPT );
	if [[ $f == *reg*-*.fa && $f != *reg*-*.fa.align.fa* && $f != *reg*-*.fa.masked && $f != *reg*-*.fa.out && $f != *reg*-*.fa.tbl && $f != *reg*-*.fa.cat* && $f != *reg*-*.fa.RM.log ]]; then
		echo "     => masking $f";
		time nohup $rm $f -pa 2 -e ncbi -a -x $cmd > $f.RM.log &
		wait ${!}
	fi
END_OF_BASH_SCRIPT
	}

# 	my @files = `ls $checks`;
# 	foreach my $f (@files) {
# 		chomp($f);
# 		$f = $checks."/".$f; #put back the path		
# 		if ($f !~ /\.align\.fa$/) { #fasta file, mask it
# 			my $cmd = "";
# 			($sp eq "na")?($cmd = "-lib $lib"):($cmd = "-sp $sp");
# 			`nohup $rm $f -pa 1 -e ncbi -a -x $cmd > $f.RM.log &`; #this will use 5 cpus per masking
# 		}
# 	}

	#Now call the script to rewrite all that in uc
	print STDERR " --- Rewriting all alignments with lower cases\n" if ($v);
	`perl $uc $checks`;
	return 1;
}








