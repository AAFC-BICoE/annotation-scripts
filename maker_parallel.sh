#!/bin/bash
#$ -S /bin/bash

# This script will parallelize maker by splitting the given fasta into individual contigs 
# This script will take advantage of SGE's task arrays
# Usage: maker_parallel.sh <file input> <working directory> '<ID field separator>' [priority] 
# 
# Parameters:
# <file input> 
#	is simply the genome/multi-fasta file you want to work with. This script will copy
# 	The fasta file into the working directory, if it is not already there. This protects
# 	source data from any clobber. The path to the file can be relative.
#
# <working directory> 
#	The folder where all the output will be found. In addition to the working directory,
# 	there will be two additional folders inside. The script will make a data folder to 
#	hold the analysis of each contig and a log folder to contain all of qsubs output. 
#	This leaves the working directory to contain the source fasta file, the folders,
#	Maker's configuration files, and a datastore index.
#
# <ID field separator> 
#	A single quoted string that separates information in the fasta file. Each sequence
#	in a fasta file has an indentification line that starts with a ">". This line contains 
#	the name and additional information. The information is usually separated by some 
#	character. However, the delimiting character isn't standard, so this needs to be given
#	to the script. PLEASE note that this should be in single quotes. (e.g. ' ')
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

# This function is used for each individual task. Rather than screating a second script, we reuse the script
qsub_maker () {
	num=$(( $SGE_TASK_ID - 1 )) # Arrays begin their index at 0, so we must shift the numbers accordingly

	cd $new_dir #Just makes sure we're in the correct directory
	
	# Bash does not allow arrays to be exported (due to safety)
	# We must retrieve our file name from a index file (See Line 118)
	file_array=( $(cat logs/index) )

	# If data exists as a file, we throw a warning to prevent clobber
	if [[ -e data && (! -d data) ]]; then
		echo "There is a file named data in the working directory. Please move it to resume and restart"
		exit 0;
	fi

	mkdir -p data
	
	cd data
	mkdir -p ${file_array[$num]} #This will allow us to make each individual folder

	# Populates/Gathers files to begin maker
	mv $new_dir/${file_array[$num]}.fasta ${file_array[$num]}/
	cp ../maker_{opts,bopts,exe}.ctl ${file_array[$num]}/
	cd ${file_array[$num]}

	# Adjusts the genome to point to the new file.
	sed -i s%^genome=.*%genome=$new_dir/data/${file_array[$num]}/${file_array[$num]}.fasta% maker_opts.ctl
	# Adjusts cpu usage as requested by user
	sed -i s%^cpus=.*%cpus=$cpu% maker_opts.ctl

	/isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker

	# Edits and copy datastore file to the global datastore file.
	while read line
	do
		name="$(echo $line | cut -d ' ' -f 1)"
		echo $name
		loca="$new_dir/data/${file_array[$num]}/${file_array[$num]}.maker.output/$(echo $line | cut -d ' ' -f 2)"
		echo $loca
		stat="$(echo $line | cut -d ' ' -f 3)"
		echo $stat
		printf "$name\t$loca\t$stat\n" >> index.tmp
	done < ${file_array[$num]}.maker.output/${file_array[$num]}_master_datastore_index.log

	cat index.tmp >> $new_dir/$filename\_master_datastore_index.log
	rm index.tmp
}

# MAIN FUNCTION BEGINE

# Checks if we are in qsub, if so begin running maker
if [[ ! -z $SGE_TASK_ID ]]; then
	qsub_maker
	exit 0;
fi

#A safe way to catch errors early. Ends the script if a command fails.
set -e

# Outputs the standard help/usage info
if [[ -z $1 || -z $2 ]]; then
        echo "Usage: $0 <fasta_file> <working_directory> '<identifier line delimiter>' [priority] [cpus]"
        echo "This script will split a given fasta file into n parts. Each file will be in their own folder. The maker configuration files must be in the working directory."
	echo "Note the single quotes around the delimiter!"
        exit 0;
fi

# Checks for valid file
if [[ ! ($1 == *.fa*) || ($1 == *.genome*) ]]; then
	echo "Unrecognized fasta file. File must have .fa* or .genome ending."
	exit 0;
fi

# Checks if you provided maker's configuration file.
if [[ ! (-a $2/maker_bopts.ctl && -a $2/maker_exe.ctl && -a $2/maker_opts.ctl) ]]; then
	echo "Cannot find maker configuration files. You can generate blank configuration files by using maker's -CTL flag."
	exit 0;
fi

# Assigns a default priority value
if [[ -z $4 ]]; then
	priority=0
else
	priority=$4
fi

# Assigns a default priority value
if [[ -z $5 ]]; then
        cpu=1
else
        cpu=$5
fi

printf "Creating initial variables... "

# Creates a copy of the genome in the working directory if it's not already there
# This is to prevent clobber and write permission issues
if [[ ! -e $2/$(basename $1) ]]; then
	cp $1 $2/
fi

# A two step process to just get the basename of a file w/o extension
filename=$(basename "$1")
filename="${filename%.*}"

# Aquires absolute paths of script and working directory
cd $(dirname $2)
export new_dir="`pwd -P`"

cd $(dirname $0)
export script_path="`pwd -P`" 

printf "done!\nNow counting contigs... "

cd $new_dir

num_contig=`grep -c "^>" $(basename $1)` #Just counts for the number of identifier lines.

printf "done!\n"

echo There are $num_contig contigs in this file. Now splitting file and creating index.

# We will use a maker script to split the genome file.
# The new names of the files will come from the identifier.
/isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/fasta_tool --split $(basename $1)

mkdir -p logs # Create a folder to store all the logs
cd logs

# As of this moment, there is no safe way to export arrays to a child process.
# As a workaround, we will dump this array to a file and read it again to rebuild.
IFS=$'\n'
export file_array=( $(grep "^>" ../subset.fasta | cut -d"$3" -f 1 | cut -d'>' -f 2) )

touch index

for i in "${file_array[@]}"
do
	echo $i >> index
done

qsub -N "Maker_Parallel" -cwd -p $priority -j y -r y -v new_dir=$new_dir,filename=$filename,cpu=$cpu -V -pe orte $cpu -t 1-$num_contig:1 $script_path/$(basename $0)

