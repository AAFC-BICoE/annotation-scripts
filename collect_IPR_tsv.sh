#!/bin/bash

#
if [ $# -eq 0 ]; then
	echo "Usage: $0 -d <MAKER datastore log file> [-o <output base> | -r]"
	echo "This script will go through the datastore directories and collate InterProScan tsv files in the current directory."
	echo "Use -r to run InterProScan instead."
	exit 0;
fi

RUN_INTERPRO=0

while [ ! -z $1 ]; do
	case $1 in
		-d) shift; log=$1; shift;;
		-o) shift; output=$1; shift;;
		-r) shift; RUN_INTERPRO=1;;
		*)  echo "Invalid Input."; echo "Usage: $0 -d <MAKER datastore log file> [-o <output base> | -r ]"; exit 1;;
	esac
done

#Checks if log can be found.
if [[ ! -e $log ]]; then
	echo "Cannot find file $log.";
	exit 1;
fi

while read line; do
	#Split line
	name=$(echo $line | cut -f1 -d ' ')
	location=$(echo $line | cut -f2 -d ' ')
	state=$(echo $line | cut -f3 -d ' ')
	
	#Only used the ones that are sucessfully completed.
	if [[ $state == "FINISHED" ]]; then
		echo "$location""$name.tsv"
		if [ -e "$location""$name.tsv" ]; then
			echo "$location""$name.tsv"
		fi

		if [ "$RUN_INTERPRO" -eq "1" ]; then
			cd $location
			file=`ls "*.all.maker.proteins.fasta" &> /dev/null`
			if [ -z $file ]; then
				echo "No protein sequences found at $name, skipping..."
			else
				qsub -N "InterProScan" -cwd -pe smp 2 -q all.q -b y /isilon/biodiversity/pipelines/interproscan-5/interproscan.sh -i $file -f HTML,GFF3,TSV -goterms -iprlookup -pa -t p
			fi
		else
			cat "$location""*.tsv" >> "${output:-output}.tsv"
		fi
	fi

done < $log
