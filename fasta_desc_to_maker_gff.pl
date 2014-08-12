#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $options = {};
GetOptions( $options,
    'fasta|f=s',
    'gff_in|i=s',
    'gff_out|o=s',);

my %name2desc = ();

sub parse_fasta
{
    my $fname = $options->{fasta};
    open(FASTA, '<', $fname) or die "Error: could not open file $fname\n";
    while (my $line = <FASTA>) {
        if ($line =~ /^>(.*)$/) {
            my @fields = split(/\|/, $1);
            if (scalar @fields == 3) {
                my $name = $fields[0];
                $name =~ s/^\s+//;
                $name =~ s/\s+$//;
                my $desc = $fields[2];
                $desc =~ s/^\s+//;
                $desc =~ s/\s+$//;
                $name2desc{$name} = $desc;
                #print "Setting $name to \"$desc\"\n";
            } else {
                print "Warning: wrong number of fields in fasta line:\n$line";
            }
        }
    }
}

sub add_gff_descriptions
{
    my $gffinfile = $options->{gff_in};
    my $gffoutfile = $options->{gff_out};
    open (GFFIN, '<', $gffinfile) or die "Error: couldn't open file $gffinfile\n";
    open (GFFOUT, '>', $gffoutfile) or die "Error: couldn't open file $gffoutfile\n";
    while (my $line = <GFFIN>) {
        if ($line =~ /Name=(BDET_[0-9]+);/) {
            my $name = $1;
            $name =~ s/^\s+//;
            $name =~ s/\s+$//;
            my $desc = ($name2desc{$name} ? $name2desc{$name} : "");
            if ($desc) {
                chomp $line;
                $line .= "Desc=\"$desc\"";
                print GFFOUT $line . "\n";
            } else {
                print STDERR "Warning: Found line with match but no description:\n$line";
                print GFFOUT $line;
            }
        } else {
            print GFFOUT $line;
        }
    }
    close (GFFIN);
    close (GFFOUT);
}

$options->{fasta} and $options->{gff_in} and $options->{gff_out} or die "Usage: $0 -f <fasta file with descriptions> -i <input gff> -o <output gff>\n";
parse_fasta
my $desc = $name2desc{"BDET_00939"};
add_gff_descriptions;
