#!/usr/bin/perl -w
#------------------------------------------------------------------------------
# Author  :  Aurelie Kapusta - https://github.com/4ureliek
# email   :  4urelie.k@gmail.com
# PURPOSE :  Util for DelGet pipeline
#------------------------------------------------------------------------------
use strict;
use warnings;
use Carp;
use Getopt::Long;

#------------------------------------------------------------------------------
#--- CHANGELOG & USAGE --------------------------------------------------------
#------------------------------------------------------------------------------
my $VERSION = "1.0";
my $SCRIPTNAME = "DelGet_cat-outputs.pl";
my $CHANGELOG = "
# UPDATES:
#	- v1.0 = 16 Mar 2018
\n";
my $USAGE = "\n USAGE [v$VERSION]:
   Typically:
     perl script.pl Del
     
   Where Del = folder where Deletions.X directories are
     
   Other options:
     -l (BOOL) prints change log / updates
     -v (BOOL) prints the pipeline version
\n";

#------------------------------------------------------------------------------
#--- PREP & MAIN --------------------------------------------------------------
#------------------------------------------------------------------------------
my ($CHLOG,$IFV);
GetOptions ('l' => \$CHLOG, 
            'v' => \$IFV);
die "\n   script $SCRIPTNAME v$VERSION\n\n" if ($IFV);
die $CHANGELOG if ($CHLOG);

my $PATH = shift or die "$USAGE $!\n" ;
$PATH =~ s/\/$//;
die "\t    ERROR - $PATH does not exist or is not a directory?\n\n" if (! -d $PATH);

# Gather the  files and store stuff to concat
#	./Del/Deletions.0/__1-30.RESULTS.tab
#	./Del/Deletions.0/__31-x.RESULTS.tab

#concatenate
concat_and_print_results("1-30");
concat_and_print_results("31-x");
exit;

#------------------------------------------------------------------------------
#--- SUBROUTINES --------------------------------------------------------------
#------------------------------------------------------------------------------
sub concat_and_print_results {
	my $type = shift;
	my @FILES = `ls $PATH/Deletions.*/__$type.RESULTS.tab` or confess "\t    ERROR - can't list files in $PATH/Deletions.* $!\n";
	my %data = ();
	foreach my $res (@FILES) {
		chomp $res;
		open(my $fh,"<",$res) or confess "\t    ERROR - could not open to read $res $!\n";
		while(defined(my $l = <$fh>)) {	
			chomp $l;
			next if ($l =~ /^#|^group/ || $l !~ /\w/);
			my @l = split(/\t/,$l);
			my $id = "$l[1]\t$l[2]\t$l[3]";
			$data{$id}{'c'}+=$l[4];
			$data{$id}{'n'}+=$l[5];
			$data{$id}{'a'}+=$l[6];
			#save outgroup info
			$data{$id}{'o'}=$l[9];
		}
		close $fh;
	}
	
	my $out = $PATH."__$type.RESULTS.all.tab";
	open(my $fh,">",$out) or confess "\t    ERROR - could not open to write file $out $!\n";
	print $fh "#Results from $SCRIPTNAME v$VERSION - ran on $PATH\n";
	print $fh "#type\tspecies\ttype\tgroup\tcount\tamount\taln_len\t%_of_aln\tdel_nt/10kb\tif_outroup(don't_use_the_numbers)\n";
	foreach my $id (keys %data) {
		my $c   = $data{$id}{'c'};
		my $nt  = $data{$id}{'n'};
		my $aln = $data{$id}{'a'};
		my $ifo = $data{$id}{'o'};
		$ifo =~ s/:can't_use_the_numbers//;
		my $gap_per = 0;
		$gap_per = $nt / $aln * 100;
		my $del_nt = $nt /  $aln * 10000;
		print $fh "$type\t$id\t$c\t$nt\t$aln\t$gap_per\t$del_nt\t$ifo\n";
	}
	close $fh;
	print STDERR "     => in $out\n";
	return 1;
}






