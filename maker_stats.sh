#!/bin/bash

# Given a directory containing some maker output merged gff/fasta files,
# calculate basic stats e.g. sequence (gene) counts.

indir=$1

# Get basic counts of number of sequences in each output fasta file.
get_counts()
{
    maker_dir=$1
    fasta_list=`find ${maker_dir} -maxdepth 1 -name "*maker*" -and -name "*.fa*" | sort`
    printf "File\tNumber of sequences\tSequence length sum\n"
    for fasta_file in $fasta_list; do
        num_seqs=`grep -c '^>' ${fasta_file}`
        seq_len=`cat ${fasta_file} | perl -ne 'BEGIN { $sum=0; } unless( /^>/ ) { chomp; $sum += length($_); } END { print "$sum\n"; }'`
        printf "${fasta_file}\t${num_seqs}\t${seq_len}\n"
    done
}
        