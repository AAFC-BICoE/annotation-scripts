#!/bin/bash
#$ -S /bin/bash

# This script will parallelize maker by splitting the given fasta into individual contigs 
# This script will take advantage of SGE's task arrays
# Usage: maker_parallel.sh -f <file input> -w <working directory> -d <ID field separator> -p [priority]
# -c [cpu]
# 
# Parameters:
# -f <file input> 
#	is simply the genome/multi-fasta file you want to work with. This script will copy
# 	The fasta file into the working directory, if it is not already there. This protects
# 	source data from any clobber. The path to the file can be relative.
#
# -w <working directory> 
#	The folder where all the output will be found. In addition to the working directory,
# 	there will be two additional folders inside. The script will make a data folder to 
#	hold the analysis of each contig and a log folder to contain all of qsubs output. 
#	This leaves the working directory to contain the source fasta file, the folders,
#	Maker's configuration files, and a datastore index.
#
# -d <ID field separator> 
#	Each sequence in a fasta file has an indentification line that starts with a ">".
#	This line contains the name and additional information. The information is usually 
#	separated by some character. However, the delimiting character isn't standard, 
# 	so this needs to be given to the script.
#
# Optional Arguments:
# Priority
#	An integer ranging from -1023 to 0. This is given to qsub to determine how urgent a task is.
# 
# CPU
#	An integer indicating how many cpus to allocate to each task.
#
# Additional Information:
# The ideology behind the script is to replace Maker's native MPI support. Rather than relying
# on Maker (which simply creates multiple processes that reads and writes the same files), we organize
# the folder so one process will work with on contig. This makes the log file clearer when we need to
# review the output. To further replace the Maker's MPI setting, a "global" datastore index is created in
# the working directory to allow users to use Maker's gff3_merge and fasta_merge.
# 

# This function is used for each individual task. Rather than screating a second script, 
# we reuse the script. This function should not be invoked directly.
# Pre: SGE variables and SGE_TASK_ID must be initialized. Runs only when the script is is called by qsub
# In addition, the current working directory should be in logs/.
# Post: Creates a folder holding all of maker's output and a fasta file of a contig from the input genome
# file. In addition, creates a master_datastore log in the working directory.
qsub_maker () {
	# Moves the current working directory up one level.
	# This is a work around to allow all of SGE's logs to remain in logs/
	cd $new_dir
	
	# Pulls the nth sequence from the fasta file.
	awk -v RS='>' -v seq=$SGE_TASK_ID 'NR>=(seq+1)&&NR<=(seq+1){print ">"$0}'  $new_dir/$(basename $from) > contig_$SGE_TASK_ID.fasta

	# Pulls the name from the fasta file
	contig_name="`head -n 1 contig_$SGE_TASK_ID.fasta`"

	if [[ -n "$delimiter" ]]; then
		contig_name=`echo ${contig_name%%$delimiter*} |  cut -d'>' -f 2`
	else
		contig_name=`echo $contig_name | cut -d'>' -f 2`
	fi

	# This decides which folder structure to use. If there's less than 300 contigs, it should be fine
	# to leave it all in one folder. Otherwise, we'll group into 50 folders.
	# We will use folder_struct to hold the folder names after data.
	if [ $SGE_TASK_LAST -ge 300 ]; then
		# Determines how many contigs in each folder
		per_folders=$(( $SGE_TASK_LAST / 50 ))
		
		# Determines which folder this contig will go into
		folder=$(( $SGE_TASK_ID / $per_folders ))

		# This is just to have the last folder matching the number of the last contig.
		if [[ ( $(( $SGE_TASK_ID % 50 )) -ne 0 ) && ( $folder -eq 50 ) ]]; then
			folder_struct="$((per_folders * 50))-$SGE_TASK_LAST/$contig_name"
		else
			min_num=$(($folder * $per_folders))
			max_num=$(( ($folder+1) * $per_folders - 1 ))
			folder_struct="$min_num-$max_num/$contig_name"
		fi
	else
		folder_struct=$contig_name
	fi
	
	echo $folder_struct

	if [[ ! -d data/$folder_struct ]]; then
		# Creates the folders, if not already there.
		mkdir -p data/$folder_struct
		# Get the fasta file to the new folder
		mv contig_$SGE_TASK_ID.fasta data/$folder_struct/$contig_name.fasta
	else
		rm contig_$SGE_TASK_ID.fasta
	fi	

	cd data
	
	if [[ ! ( -e $folder_struct/maker_opts.ctl && -e $folder_struct/maker_bopts.ctl && -e $folder_struct/maker_exe.ctl ) ]]
	then
		# Gathers Maker config files
		cp ../maker_{opts,bopts,exe}.ctl $folder_struct/
	
		cd $folder_struct
		# Adjusts the genome to point to the new file.
		sed -i s%^genome=.*%genome=$new_dir/data/$folder_struct/$contig_name.fasta% maker_opts.ctl
		# Adjusts cpu usage as requested by user
		sed -i s%^cpus=.*%cpus=$cpu% maker_opts.ctl
	
		# Changes temp directory to avoid NFS errors. (Bug 2183 Comment 19)
		sed -i s%^TMP=.*%TMP=/state/partition1/% maker_opts.ctl
	else
		cd $folder_struct
	fi

	maker # Run Maker. It should have been installed so that it gets added to the PATH.

	# Edits and copy datastore file to the global datastore file.
	while read line
	do
		name="$(echo $line | cut -d ' ' -f 1)"
		loca="$new_dir/data/$folder_struct/$contig_name.maker.output/$(echo $line | cut -d ' ' -f 2)"
		stat="$(echo $line | cut -d ' ' -f 3)"
		printf "$name\t$loca\t$stat\n" >> index.tmp
	done < "$new_dir/data/$folder_struct/$contig_name.maker.output/$contig_name""_master_datastore_index.log"

	cat index.tmp >> "$new_dir/$filename""_master_datastore_index.log"
	rm index.tmp
	
	# Optionally Run InterProScan
	if [[ $stat = "FINISHED" && ! -z $interpro && ! -e $loca/$name.tsv ]]; then
		cd $loca
		$interpro -b $name -f TSV,XML,GFF3,HTML -goterms -iprlookup -pa -i $name.maker.proteins.fasta #run interpro
		ipr_update_gff $name.gff $name.tsv #Merge results back to gff file
	fi

	# Explicitly define UTR regions and start and stop codons
	add_utr_start_stop_gff $name.gff tmp-$name.gff
	mv tmp-$name.gff $name.gff

	cd $new_dir

	# Collates all the results into a genome level file.
	gff3_merge -d "$filename""_master_datastore_index.log" -g -o $filename.genes_only.gff3
	gff3_merge -d "$filename""_master_datastore_index.log" -o $filename.maker.all.gff3
	fasta_merge -d "$filename""_master_datastore_index.log"
}

# MAIN FUNCTION BEGIN

# Checks if we are in qsub, if so begin running maker
if [[ ! -z $SGE_TASK_ID ]]; then
	qsub_maker
	exit 0;
fi

#A safe way to catch errors early. Ends the script if a command fails.
set -e

usage="\nUsage: $0 -f <fasta_file> -w <working_directory> [-i <interpro location>-d <identifier line delimiter> -p [priority] -c [cpus]] \n\nThis script will split a given fasta file so each sequence is in their own fasta file. Each file will be in their own folder. \nThe maker configuration files must be in the working directory. \nYou may specify to run a interproscan by using -i <location>. \n\nFor more information, read the documentation at the top of this script."

if [ $# -eq 0 ]; then
        echo -e $usage
        exit 0;
fi

while [ ! -z $1 ];
do
	case $1 in
		-f) shift; file=$1; shift;;
		-w) shift; new_dir=$1; shift;;
		-i) shift; interpro=$1; shift;;
		-d) shift; delimiter="$1"; shift;;
		-p) shift; priority=$1; shift;;
		-c) shift; cpu=$1; shift;;
		*) echo -e $usage; exit 0;;
	esac
done

# Checks for valid file
if [[ ! ($file == *.fa*) || ($file == *.genome*) ]]; then
	echo "Unrecognized fasta file. File must have .fa* or .genome ending."
	exit 1;
fi

# Checks if you provided maker's configuration file.
if [[ -z $new_dir || ! (-e $new_dir/maker_bopts.ctl && -e $new_dir/maker_exe.ctl && -e $new_dir/maker_opts.ctl) ]]; then
	echo "Cannot find maker configuration files. You can generate blank configuration files by using maker's -CTL flag."
	exit 1;
fi

if [[ -n $interpro && ! -f $interpro ]]; then
	echo "Cannot find $interpro. Please ensure that the path is valid."
	exit 1;
fi

# Assigns a default priority value
if [[ -z $priority ]]; then
	priority=0
fi

# Assigns a default priority value
if [[ -z $cpu ]]; then
        cpu=1
fi

printf "Creating initial variables... "

# Creates a copy of the genome in the working directory if it's not already there
# This is to prevent clobber and write permission issues
if [[ ! -e $new_dir/$(basename $file) ]]; then
	cp $file $new_dir/
fi

# A two step process to just get the basename of a file w/o extension
filename=$(basename "$file")
filename="${filename%.*}"

# Aquires absolute paths of script and working directory
cd $(dirname $new_dir)
export new_dir="`pwd -P`"

cd $(dirname $0)
export script_path="`pwd -P`" 

printf "done!\nNow counting contigs... "

cd $new_dir

num_contig=`grep -c "^>" $(basename $file)` #Just counts for the number of identifier lines.

printf "done!\n"

echo There are $num_contig contigs in this file.

if [[ ( -d logs ) && -n "`ls logs`" ]]; then
	set +e
        mkdir -p old_logs
        mv logs/* old_logs/
	mv $filename"_master_datastore_index.log" old_logs/
	set -e
else
	mkdir -p logs # Create a folder to store all the logs
fi

cd logs

qsub -N "Maker_Parallel" -V -pe smp $cpu -t 1-$num_contig:1 -tc 10 -cwd -q all.q -p $priority -b y -v new_dir=$new_dir,filename=$filename,cpu=$cpu,from=$file,delimiter="$delimiter",interpro=$interpro $script_path/$(basename $0)

