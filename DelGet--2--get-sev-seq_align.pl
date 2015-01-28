#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie K
# version :  3.2 (see below)
# email   :  4urelie.k@gmail.com  
# Pupose  :  Specifically written to extract orthologous sequences from several genomes + align them
#				input file needs to be as: ID \t scaffold/chr \t start \t end \t strand \t SpeciesID \t genome_location
#				Typically:
#
#					reg1     scaffold/chr     start     end     +     Mluc     /data/genomes_to_align/Myotis_lucifugus/Myotis_7x.fasta
#					reg1     scaffold/chr     start     end     +     Mdav     /data/genomes_to_align/Myotis_davidii/fasta/Mdavidii.assembly2012.fa
#					reg1     scaffold/chr     start     end     -     Efus     /data/genomes_to_align/Eptesicus_fuscus/Eptesicus_fuscus_assembly1.0.fasta
#					reg2     scaffold/chr     start     end     +     Mluc     /data/genomes_to_align/Myotis_lucifugus/Myotis_7x.fasta
#					reg2     scaffold/chr     start     end     -     Mdav     /data/genomes_to_align/Myotis_davidii/fasta/Mdavidii.assembly2012.fa
#					reg2     scaffold/chr     start     end     +     Efus     /data/genomes_to_align/Eptesicus_fuscus/Eptesicus_fuscus_assembly1.0.fasta
#
#				Sequences of region1 will be aligned
#				Note that if strand is -, reverse complement will be extracted
#######################################################
# UPDATES
#	- v1.0 = Apr 2013
#   - v2.0 = May 2013
#		- check step to avoid extracting (and therefore aligning) if region is too big
#		- do not create new output folder anymore [to allow checking for file existence]
#		- check step to avoid aligning is muscle.fa file exists
#   - v2.1 = Jul 2013
#		- change of $multip variable value (dep. if 10kb or 25kb)
#		- modification on fasta seq names match to extract sequences [spe primates or carnivoras]
#		- even if no extraction, do alignment
#   - v3.0 = 14 Oct 2013
#		- reintegration in pipeline => config file is used
#		- modification of extraction step (dealing with chr names)
#   - v3.1 = 17 Oct 2013
#		- muscle parameter change for very long sequences, to avoid error *** ERROR ***  MSA::GetLetter
#   - v3.2 = 18 Oct 2013
#		- Kalign used instead of muscle for large sequences [muscle still not working]
#######################################################
use warnings;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	unshift(@INC, "$BIN/Lib");
}
use Array::Unique;
my $version = "3.2";

my $usage = "\nUSAGE:
	perl <scriptname.pl> <config_file> <path> <inputfile>
	
	This script is part of a pipeline, type \"perl DelGet.pl -help\" for more details about it + check config file
	
	path corresponds to the /path_to/Deletions.X folder where previous outputs were generated
	inputfile is as: ID \t scaffold/chr \t start \t end \t strand \t speciesID \t genome_location
	
	Script will extract sequences corresponding to the same region in different genomes + align them using muscle [ie align orthologous regions]\n\n";

my $ArgN = @ARGV;
if ($ArgN != 3) {
	die $usage;
}
my $config_file = shift or die "$usage";
my $pathtemp = shift or die "$usage";
my $file = shift or die "$usage";

#################################################################
# Variables and names
#################################################################
# Initialize needed configuration file variables [see config_file]
our $anch_dist; 
our $multip;
our $ALNSOFT; #will be muscle or Kalign, based on $anch_dist

# Load Configuration file
do "$config_file";

my $max = $multip*$anch_dist;

#LOG FILE
my $log = "$pathtemp/_DelGet--2--get-sev-seq_align.log";
open(LOG, ">$log") or die "\t    ERROR - can not create log file $log $!\n";

# check for variables, kill this script and pipeline if some are missing
die print LOG "\t!! Some variables are not defined in config file\n" if ((! $MUSCLESOFT) || (! $anch_dist) || (! $multip));

# if OK, then carry on
print LOG "\n--- Script (v$version) started with following muscle location (can be changed inside Config file):
    Muscle software = $MUSCLESOFT
    Max sequence length = $anch_dist x $multip = $max nt\n\n";

############################################################
# Create directory to work in
############################################################
# Create output directory
my $path = "$pathtemp/_ExtractAlign";
mkdir ($path, 0755) or die print LOG "\t !! Couldn't mkdir $path $!" unless (-e $path);

####################################################################
# EXTRACT SEQUENCES
####################################################################
my $i = 0;
print LOG "\n--- extracting sequences...\n";
my @RegList = ();
tie @RegList, 'Array::Unique';
my @NoAlign = ();
tie @NoAlign, 'Array::Unique';
my %regcheck = ();
open(POSI,"<$file") or die "\t!! can't open file $file $!\n";
LINE: while(<POSI>){
	next LINE if (($_ !~ /\w/) || ($_ =~ /^#/));	#avoid blank lines and header lines
	chomp (my $line = $_);	
	my ($region,$chr,$start,$end,$strand,$species,$genome) = split(/\t/,$line);
	unless (($region) && ($chr) && ($start) && ($end) && ($strand) && ($species) && ($genome)) {
		print LOG "\t !! $line does not contain all necessary columns -> next $line\n";
		next LINE;
	}
	$strand="RevCom" if ($strand eq "-");
	#0			1	2		3	4	5		6
	#region1	chr	start	end	+	Mluc	/mnt/disk3/genomes/Myotis_lucifugus/Myotis_7x.fasta
	
	####################################################################
	# DEAL WITH GENOME
	####################################################################
	my $reindex;
	my $indexfile = "$genome.index";
	if (-e $indexfile) {
		$reindex = 0;
		#print "\t--- Genome previously indexed - Skipping indexing... (if you want to reindex the genome, delete $indexfile)\n";
	} else {
		$reindex = 1;
		print LOG "\tGenome not indexed - Indexing $genome...\n";
	}
	my $db = Bio::DB::Fasta->new($genome, -reindex=>$reindex) or die "\t !! Failed to create Bio::DB::Fasta object from $genome $!\n";
	my @dbIDs = $db->get_all_ids();
	####################################################################
	
	my $regionseq = "$path/$region.fa";
	# CHECK EXISTENCE
	if ((-f $regionseq) && (! $regcheck{$region})) {#ie if file exists despite fact this is the first extraction for this region
		print LOG "\t!! $regionseq exists => no extraction, but will still be aligned\n";
		$i = 1;
		push(@RegList,$regionseq);
		next LINE; 
	}
	$regcheck{$region}=1;
	# IF OK THEN GO AHEAD
	my $out = Bio::SeqIO->newFh(-format => 'Fasta', -file=>">>$regionseq") or die print LOG "\t !! Failed to create SeqIO FH object from $regionseq $!\n";
	
	# now deal with the | in scaffold names, problem to match
	my $expr = $chr;
	$expr =~ s/\|/\\\|/g;
	
	#name of the sequence extracted will be only gb when there is gi and gb in name (otherwise too long)
	my $chrmod = $chr;
	$chrmod =~ s/^gi\|.*\|gb/gb/;
	
	# prepare new name
	my $newId = ($species."_".$chrmod.":".$start."-".$end."_".$strand);
	
	#now extract
	if ("@dbIDs" =~ m/(\S*)($expr)(\S*)/){
		my $subSeq = $db->seq($chr,$start,$end);	# extract target sequence
		my $seqobj = Bio::Seq->new( -display_id => $newId, -seq => $subSeq); # create object with target
		my $len = $seqobj->length;
		$seqobj = $seqobj->revcom if ($strand eq "RevCom"); # reverse complement if hit was antisense
		print LOG "\t\t$region - $species => $chr:$start-$end => $len\n";
		if ($len >= $max){
			print LOG "\t!! seq for $species not extracted, its length > max ($len > $max)\n\t\t=> $region.fa won't be aligned\n";
			push(@NoAlign,$regionseq);
			$i = 1;
		} else {
			print $out $seqobj;	# print it out (in fasta format)
			push(@RegList,$regionseq);
		}	
	}
	else {
		print LOG "\t!! $chr was not found in $genome\n";
		$i = 1;
	}
}
close POSI;
print LOG "--- done\n";

####################################################################
# ALIGNING SEQUENCES
####################################################################
print LOG "\n--- Aligning sequences...\n";

# print LOG "    From this list:\n";
# foreach my $al (@RegList) {
	# print LOG "\t$al\n";
# }
# print LOG "\n";

print LOG "    Except for the following ones:\n" if ($#NoAlign >= 1);
foreach my $noal (@NoAlign) {
	print LOG "\t$noal\n";
}
print LOG "\n";

foreach my $seqtoalign (@RegList) {
	my $alnout = "$seqtoalign.align.fa";
	unless ((-f $alnout) || ("@NoAlign" =~ /$seqtoalign/)){
		print LOG "    Aligning sequences in $seqtoalign...\n";
		my $alignlog = "$seqtoalign.align.log";
		open ALN, ">$alnout" or die print LOG "\t!! can not create alnout file $alnout $!\n";
		system "$ALNSOFT -in $seqtoalign -out $alnout -log $alignlog -quiet -verbose" || die print LOG "\t!! Failed running muscle: $!\n" if ($ALNSOFT =~ /[Mm]uscle/);
		#system "$ALNSOFT -in $seqtoalign -out $alnout -log $alignlog -quiet -verbose -maxiters 2 -diags1 -sv" || die print LOG "\t!! Failed running muscle: $!\n" if ($anch_dist > 40000);
		system "$ALNSOFT -gpo 80 -gpe 3 -tgpe 3 -bonus 0 -format fasta -quiet -in $seqtoalign -out $alnout" || die print LOG "\t!! Failed running kalign: $!\n" if ($ALNSOFT =~ /[Kk]align/);
			#defaults are: gap open penalty      -gpo   = 217.00000000
			#			   gap extension penalty -gpe   =  39.40000153
			#			   terminal gap penalty  -gpe   = 292.60000610
			#			   bonus                 -bonus = 283.00000000
			
		close ALN;
		unlink $alignlog;
	} else {
		$i = 1;
		print LOG "\t!! $alnout exists => Not realigning $seqtoalign\n" if (-f $alnout);
		print LOG "\t!! $seqtoalign had non extracted sequences => Not aligning\n" if ("@NoAlign" =~ m/(\S*)($seqtoalign)(\S*)/);
	}
}
print LOG "--- done\n";
####################################################################
if ($i == 1) {
	print "\t    !! There were some errors, please check $log\n";
}
print LOG "\n--- Finished, see files in $path\n\n";
exit;

## NB for memory:
#sub get_all_ids {
# grep {!/^__/} keys %{shift->{offsets}}
#}



