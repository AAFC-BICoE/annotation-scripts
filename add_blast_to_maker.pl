use strict;
use warnings;

my $usage = "
add_blast_to_maker.pl maker.gff ncbi-blast.outfmt6.txt blastdb {nucl,prot} threshold_value new_file.gff

This script will take a blast report and attempt to create new GFF3 features from 
matches of at least the given value. Prints to STDOUT. Choose nucl or prot to specify the
database type!

The threshold value should be an integer (like 70) to represent 70%. We will query the 
blastdb using blastdbcmd for the description and list any additional names. The blast 
results must be in the following format:

qseqid sseqid pident evalue bitscore qstart qend

You can create a blast report with more than these columns for your own use. You can then generate 
your blast output and use the following awk one-liner to generate the output for this script.

awk \'BEGIN {FS=\"\\t\"}; {print \$1,FS,\$2,FS,\$3,FS,\$11,FS,\$12,FS,\$7,FS,\$8,FS,\$13}\' blast_results.txt

You can change the numbers to match the expected column number of your output. 

";

my ( $gff_file, $blast_file, $blastdb, $dbtype, $threshold, $output ) = @ARGV;

# Ensure that the correct parameters are given before proceding
die $usage unless $blastdb && (-e $blast_file ) && (-e $gff_file ) && $threshold && $dbtype && $output;

if ( -e $output ) {
	die "$output already exists. Will not clobber.\n";
}

#----------------------------------------
#Parse blast_file

#my $exe = which('blastdbcmd') or die "Can't find blastdbcmd\n";
my $exe = `which blastdbcmd` or die "Can't find blastdbcmd\n";
chomp $exe;

open(my $IN, '<', $blast_file) or die "Can't open $blast_file for reading\n";

print "Now parsing the blast report...\r";

#Calculate line numer for progress tracking
my $num = (split(" ", `wc -l $blast_file`))[0];
my $count = 0;

# %rna is a hash of arrays that maps
# $qseqid => [ \@hits ]
#
# @hits is an array of arrays of at least size 1.
# @hits = ( [ \@entry ] )
#
# @entry is an array/tuple of mixed entries
# @entry = ([ @gi_list ], [ @description_list ], $pident, $evalue, $bitscore, $qstart, $qend)
# 
# @database_ref and @description_list are both arrays of strings.
my (%rna, @hits, @entry);

while (<$IN>) {
	$count = $count + 1;
	#Print progress every so often
	if ( $count % 500 == 0 ) {
		no integer;
		my $progress = sprintf("%.0f", ($count/$num *100));
		print "Now parsing the blast report...$progress%\r";
	}

	chomp;	
	my ($qseqid, $sseqid, $pident, $evalue, $bitscore, $qstart, $qend) = split(/\t/);
	
	use integer;
	#If the blast hit is less than the threshold value, ignore
	if ( $pident <= $threshold ) {
		next;
	}

	#Don't assume that the file is grouped by seqid. So, we'll use a hash to keep track of the entries.
	if ( exists $rna{$qseqid} ) {
		@hits = @{ $rna{$qseqid} };
	} else {
		@hits = ();
	}

	#Find some form of ID to look in database
	my $id = (split(/\|/, $sseqid))[1];
	#Get fasta identifiers (This is probably the slowest part)
	#Should look like this: ref|WP_003131952.1|#!#30S ribosomal protein S18 [Lactococcus lactis]
	my @fasta_identifiers = split("\n", `$exe -db $blastdb -dbtype $dbtype -entry $id -outfmt \"%i\#\!\#%t\"`);
	
	#Clear the arrays for new entries
	my @database_ref = ();
	my @description_list = ();

	#For each identifier, we will create a new feature.
	#For now, keep track of the id and descriptions, worry about formatting later
	foreach (@fasta_identifiers) {
		push(@database_ref, (split(/\#\!\#/, $_))[0]);
		push(@description_list, (split(/\#\!\#/, $_))[1]);
	}
	
	#Add entry, update hits in the hash table
	@entry = ([ @database_ref ], [ @description_list ], $pident, $evalue, $bitscore, $qstart, $qend);
	push (@hits,  [@entry]);
	$rna{$qseqid} = [@hits];
}

print "Now parsing the blast report...Done!\n";
close $IN;

#---------------------------------------------------------------
# Read GFF3 file and add new entries as required


open( $IN, "<", $gff_file) or die "Can't open $gff_file.\n";
open( my $OUT, ">", $output ) or die "Can't write to $output.\n";

print "Now writing new GFF file...\r";

#Variables for progress tracking
$num = (split(" ", `wc -l $gff_file`))[0];
$count = 0;

while(<$IN>) {
        $count = $count + 1;
	#print progress every so often
	if ( $count % 500 == 0 ) {
	        no integer;
        	$count = $count + 1;
	        my $progress = sprintf("%.0f", ($count/$num *100));
	        print "Now writing new GFF file...$progress%\r";
	}

	#Get variables
	my ($seq, $source, $type, $start, $end, $score, $strand, $phase, $attribute) = split(/\t/);

	#print line first
	print $OUT "$_";

	if ($_ =~ /^\s*\#/ || ! $attribute || ! ( $type =~ /mrna/i || $type =~ /transcript/i) ) {
		next;
	}

	my ($id) = $attribute =~ /ID\s*=\s*(.*?);/;

	if ( $rna{$id} ) {
		@hits = @{ $rna{$id} };
	} else {
		next;
	}


	for my $entry ( 0 .. $#hits ) {
		#Assume that the size of @database_ref and @description_list is the same
		for my $i ( 0 .. $#{ $hits[$entry][0] } ) {

			#Determine which strand is the match on
			if ( $hits[$entry][5] < $hits[$entry][6] ) {
				$strand = '+';
			} else {
				$strand = '-';
			}

			#Parse description and ID information
			my ($name, @db_id);
			my $dbref = "";

			#UniProt: S-adenosylmethionine synthase OS=Ascobolus immersus PE=3 SV=1
			if ( $hits[$entry][1][$i] =~ /(.*?)\s+OS=(.*?)\s+(GN=(.*?)\s+)?PE=/ ) {
				$name = "$4 " if $4;
				$name .= "$1 " if $1;
				$name .= "($2) " if $2;
				$name .= "$hits[$entry][2]% identity, bitscore$hits[$entry][4]";
			#Genbank uniprot: RecName: Full=30S ribosomal protein S18
			} elsif ( $hits[$entry][1][$i] =~ /RecName: [A-Za-z]+=(.+?);*.*/) {
				$name = $1;
				$name .= " $hits[$entry][2]% identity, bitscore$hits[$entry][4]";
			#Genbank  normal: 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris UC509.9]
			} else {
				$name = $hits[$entry][1][$i];
				$name .= " $hits[$entry][2]% identity, bitscore$hits[$entry][4]";
			}

			#Parse the id's 
			@db_id = split(/\|/, $hits[$entry][0][$i]);
			#ref|YP_007509482.1|
			if ( $db_id[0] =~ /ref/ ) {
				$dbref = "RefSeq:$db_id[1]";
			#sp|P50304|METK_ASCIM
			} elsif ( $db_id[0] =~ /sp/ ) {
				$dbref = "Swiss-Prot:$db_id[1]";
			#could be any of the ones below
			} elsif ( $dbtype eq "prot" ) {
				$dbref = "protein_id:$db_id[1]";
			#gb|AAK06287.1|AE006448_5|
			} elsif ( ($db_id[0] =~ /gb/) && $dbtype eq "nucl" ) {
				$dbref = "GB:$db_id[1]";
			#emb|CAL99037.1|
			} elsif (  ($db_id[0] =~ /emb/) && $dbtype eq "nucl" ) {
				$dbref = "EMBL:$db_id[1]";
			#dbj|AB000282.1|
			} elsif (  ($db_id[0] =~ /dbj/) && $dbtype eq "nucl" ) {
				$dbref = "DDBK:$db_id[1]";
			}

			#Attribute example:
                        #ID=69326-augustus-gene-0.103-mRNA-1;Parent=69326-augustus-gene-0.103;Name=69326-augustus-gene-0.103-mRNA-1;Dbxref=Gene3D:G3DSA:1.10.20.10,InterPro:IPR003958,InterPro:IPR00907;
			my $attrib = join(";", ("ID=$id\_blast_match_$entry\_redundant-entry_$i", "Name=$name", "Dbxref=$dbref;"));
			
			my $feature;
			if ( $dbtype eq "prot" ) {
				$feature = join("\t", ($seq, "blast", "protein_match", $start, $end, $hits[$entry][3], $strand, ".", $attrib));
			} else {
				$feature = join("\t", ($seq, "blast", "nucleotide_match", $start, $end, $hits[$entry][3], $strand, ".", $attrib));
			}

			print $OUT "$feature\n";
		}
	}
}
print "Now writing new GFF file...Done!\n";
