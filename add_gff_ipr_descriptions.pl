#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# Parse interpro gff, pull out the description for each, then add this to
# the maker gff (maker gff must have already had ipr_update_gff run on it) 
# wherever we see the same interpro id.

my $options = {};
GetOptions($options,
    'maker_gff|m=s',
    'interpro_tsv|i=s',
    'output_gff|o=s',);

$options->{interpro_tsv} and $options->{maker_gff} and $options->{output_gff} or die "Usage: $0 -m <maker gff> -i <interpro tsv> -o <output gff>\n";

my %ipr2desc = ();

sub parse_interpro
{
    my $fname = $options->{interpro_tsv};
    open (FTSV, '<', $fname) or die "Error: couldn't open file $fname\n";
    while (my $line = <FTSV>) {
        chomp $line;
        my @fields = split (/\t/, $line);
        my $ipr_id = "";
        my $desc = "";
        if (scalar @fields >= 13) {
            $ipr_id = $fields[11];
            $desc = $fields[12];
            
            $ipr_id =~ s/^\s+//;
            $ipr_id =~ s/\s+$//;
            $ipr2desc{$ipr_id} = $ipr_id . "_desc=\"" . $desc . "\";";
            #print $ipr2desc{$ipr_id} . "\n";
        }
    }
    close (FTSV);
}

sub update_gff
{
    my $fname_in = $options->{maker_gff}; 
    my $fname_out = $options->{output_gff};
    my %missing = ();
    open (INGFF, '<', $fname_in) or die "Error: couldn't open file ${fname_in}\n";
    open (OUTGFF, '>', $fname_out) or die "Error: couldn't open file ${fname_out}\n";
    while (my $line = <INGFF>) {
        chomp $line;
        my $desc_str = "";
        foreach ($line =~ /(IPR[0-9]+)/g) {
            my $desc = $ipr2desc{$1};
            if ($desc) {
                $desc_str .= $desc;
            } else {
                $missing{$1}++;
            }
        }
        $line .= $desc_str;
        print OUTGFF $line . "\n";
    }
    if (scalar keys %missing > 0) {
        print "IPRs missing descriptions:\n";
        print join ("\n", (sort keys %missing)) . "\n";
    }
    close (INGFF);
    close (OUTGFF);
}

parse_interpro;
update_gff;

