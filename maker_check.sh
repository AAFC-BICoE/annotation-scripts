#!/bin/bash

## ./maker_check.sh <working directory> 
## This script will check the input files in MAKER control files to ensure that all paths are valid
## Exit codes: 	64 -> Invalid File
## 		65 -> Invalid string

. $1/maker_opts.ctl
. $1/maker_bopts.ctl
. $1/maker_exe.ctl

if [ ! -s "$genome" ]; then
	echo "Genome file is empty or does not exist: $genome"
	exit 64;
fi

if [[ ! ( "$organism_type" != "eukaryotic" || "$organism_type" != "prokaryotic" ) ]]; then
	echo "Invalid organism type. Given $organism_type"
	exit 65;
fi

if [[ "$genome_gff" && -s "$genome_gff" ]]; then
	if [[ ! ( ( $est_pass -eq 0 || $est_pass -eq 1 ) && ( $altest_pass -eq 0 || $altest_pass -eq 1 ) && ( $protein_pass -eq 0 || $protein_pass -eq 1 ) && ( $rm_pass -eq 0 || $rm_pass -eq 1 ) && ( $model_pass -eq 0 || $model_pass -eq 1 ) && ( $pred_pass -eq 0 || $pred_pass -eq 1 ) && ( $other_pass -eq 0 || $pred_pass -eq 1 ) ) ]]; then
		echo "Invalid input in Re&&nnotation settings"
		exit 65;
	fi
elif [[ "$genome_gff" && ! -s "$genome_gff" ]]; then
	echo "$genome_gff cannot be found or is empty"
	echo 64;
fi

for name in "$est" "$altest" "$est_gff" "$altest_gff" "$protein" "$protein_gff" "$rmlib" "$repeat_protein" "$rm_gff" "$snaphmm" "$gmhmm" "$fgenesh_par_file" "$pred_gff" "$model_gff" "$other_gff"; do
	if [[ "$name" && ! -s "$name" ]]; then
		echo "$name cannot be found or is empty"
		exit 64;
	fi
done

for name in "prok_rm" "est2genome" "protein2genome" "unmask" "alt_splice" "always_complete" "map_forward" "keep_preds" "single_exon" "clean_try" "clean_up"; do
	if [[ ! ( "$( echo $name )" -eq 0 || "$( echo $name )" -eq 1 ) ]]; then
		echo "$name has an invalid input. Received $(echo $name)."
		exit 65;
	fi
done

if [ "$TMP" != "/state/partition1/" ]; then
	echo "Bad TMP directory. Should be /state/partition1/";
	exit 65;
fi
