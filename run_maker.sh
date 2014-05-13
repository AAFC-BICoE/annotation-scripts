#!/bin/bash

usage() { echo "Usage: $0 [-s <start contig size> -e <end contig size> -f <function>" 1>&2; exit 1; }

cstart=
cend=
func=

while getopts "s:e:f:" opt; do
    case "${opt}" in
        s)
            cstart=${OPTARG}
            ;;
        e)
            cend=${OPTARG}
            ;;
        f)
            func=${OPTARG}
            ;;
    esac
done

if [[ $cstart && $cend && $func ]]; then
    contig_range="$cstart-$cend"
else
    usage
fi

# export contig_range="x-y"
genome_in=/isilon/biodiversity/projects/PRI_Se/assembly/clc-sendo.fasta
genome_out=Se_PRI2_assembly.fa
RNA_in=/isilon/biodiversity/users/cullisj/bug3124/CFIA_contigs.fa
RNA=CFIA_RNA_contigs.fa
proteins_in=/isilon/biodiversity/users/cullisj/bug3132/B_dendrobatidis/Batde5_best_proteins.fasta
proteins=B_dendro_JGI_proteins.fa
echo "Using contig range $contig_range"

maker1_dir=maker1
augustus1_dir=augustus1
snap_dir=snap
genemark_es_dir=genemark_es
genemark_sn_dir=genemark_sn
gm_hmm=

gensub=Se_PRI2_assembly_${contig_range}
gensub_fa=$gensub.fa
#maker_out=$maker1_dir/$gensub.maker.output
#species=synchytrium_endobioticum
augustus_species=synchytrium_endobioticum_${contig_range}

maker_bopts=maker_bopts.ctl
maker_exe=maker_exe.ctl
maker_opts=maker_opts.ctl

dir_setup()
{
    ln -s $genome_in $genome_out
    ln -s $RNA_in $RNA
    ln -s $proteins_in $proteins
}

sample_fasta()
{
    [ -e ./fastaSizes.pl ] || svn export http://biodiversity/svn/source/misc_scripts/fastaSizes.pl
    perl ./fastaSizes.pl -f $genome_out -r $contig_range -o $gensub_fa
}

check_organism()
{
    x=
    # Check that we get some hits to the target in our subset
}

replace_vars()
{
    fname=$1
    key=$2
    value=$3
    sed -i "s/^$key=[^\s\\]*\s/$key=$value #/" $fname
}

run_maker1()
{
    # Run first pass of maker
    mkdir $maker1_dir
    cp $gensub_fa $RNA $proteins $maker1_dir/
    cd $maker1_dir
    maker -CTL
    replace_vars $maker_opts genome $gensub_fa
    replace_vars $maker_opts est $RNA
    replace_vars $maker_opts protein $proteins
    replace_vars $maker_opts est2genome 1
    #replace_vars $maker_opts protein2genome 1 (only if prokaryotic)
    replace_vars $maker_opts alt_splice 1
    replace_vars $maker_opts TMP \\/state\\/partition1
    replace_vars $maker_bopts blast_type ncbi
    replace_vars $maker_exe RepeatMasker \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/RepeatMasker-open-4-0-3\\/RepeatMasker
    replace_vars $maker_exe exonerate \\/opt\\/bio\\/exonerate\\/bin\\/exonerate
    
    maker_out=$gensub.maker.output
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/gff3_merge -d $maker_out/*_master_datastore_index.log
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/fasta_merge -d $maker_out/*_master_datastore_index.log
    sed '/FASTA/q' $gensub.all.gff | sed '$d' >$gensub.all.nofa.gff
    
    # Line below creates genome.ann, genome.dna files, required by augustus, snap, ..
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker2zff -d $maker_out/*_master_datastore_index.log
    # mv $gensub.all* genome.ann genome.dna $maker_out/
    cd ..
}

setup_augustus()
{
    # Only run this if you've never set up before
    mkdir ~/augustus
    cp -r /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/config ~/augustus/
    echo "export AUGUSTUS_CONFIG_PATH=~/augustus/config" >> ~/.bashrc
    . ~/.bashrc
}

# Augustus is only for eukaryotes
# Depends on maker1 output
train_augustus1()
{
    mkdir -p augustus1
    cp $maker1_dir/$gensub.all.nofa.gff $maker1_dir/genome.dna augustus1/
    cd augustus1
    /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/scripts/gff2gbSmallDNA.pl $gensub.all.nofa.gff genome.dna 1000 genes.gb
    #above breaks if fasta is in gff file.
    /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/scripts/new_species.pl --species=$augustus_species

    # Get count of total number of annotations in genes.gb
    total_ann=`grep -c '^ORIGIN' genes.gb`
    # Take 30% for testing. This is standard for ML but not sure for annotation.
    num_train=`echo "scale=0; $total_ann * 3/10" | bc`;
    /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/scripts/randomSplit.pl genes.gb $num_train
    /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/bin/etraining --species=$augustus_species genes.gb.train  | tee etraining_initial

	/isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/bin/augustus --species=$augustus_species genes.gb.test | tee accuracy_initial
	# next line http://stackoverflow.com/questions/592620/how-to-check-if-a-program-exists-from-a-bash-script
	command -v optimize_augustus.pl >/dev/null 2>&1 || export PATH="/isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/bin:$PATH"
	/isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/scripts/optimize_augustus.pl --species=$augustus_species genes.gb.train | tee optimize_log
	/isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/bin/etraining --species=$augustus_species genes.gb.train | tee etraining_final
	/isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/bin/augustus --species=$augustus_species genes.gb.test | tee accuracy_final
}

# Without running CEGMA first here
# Depends on maker1 output
train_snap()
{
    mkdir -p $snap_dir
    cp $maker1_dir/genome.dna $maker1_dir/genome.ann $snap_dir/
    cd $snap_dir
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -gene-stats
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -validate
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom export.ann export.dna -export 1000 -plus
    mkdir params
	cd params
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/forge ../export.ann ../export.dna
	cd ..
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/hmm-assembler.pl $gensub params > $gensub.hmm
	#Note: At the last step, you can replace my-genome with a name and my-genome.hmm
	#with the file name.
}

# Train genemark for eukaryotes
train_genemark_es()
{
    mkdir $genemark_es_dir
    cp $gensub_fa $genemark_es_dir/
    #cp $maker1_dir/genome.ann $maker1_dir/genome.dna $genemark_es_dir/
    cd $genemark_es_dir
    
    /isilon/biodiversity/pipelines/maker-2.10/gene-mark-es-2.3e/gmes/gm_es.pl $gensub_fa
	# Note: If GeneMark fails, there might be something wrong with your genome.
	# If your contigs are short, try adding --min_contig 10,000 and --max_nnn 5000
    # When completed, the training file is ./mod/es.mod. Note that this is a symlink
    # to another file. If you want to move it, you should just copy the actual file.
    gm_hmm=`readlink $genemark_es_dir/mod/es.mod`
}

# Train genemark for prokaryotes
train_genemark_sn()
{
    mkdir $genemark_sn_dir
    cp $gensub_fa $genemark_sn_dir/
    cd $genemark_sn_dir
    name=gmsn
    gmsn.pl --combine --species $augustus_species -gm --name $name $gensub_fa
    # When completed, the training file is <name>_hmm_combined.mod.
}

run_maker2()
{
    mkdir $maker2_dir
    cp $gensub_fa $RNA $proteins $maker2_dir
    # cp $snap_dir/path/to/snaphmm snap.hmm
    # if [ type=eukaryote ]; then ...
    cp $gm_hmm $maker2_dir/gm_es.mod $maker2_dir/ # ln -s to be able to see target here?
    
    cp $maker1_dir/$gensub.all.gff $maker2_dir
    
    cp $maker1_dir/*ctl $maker2_dir/
    cd $maker2_dir
    #maker -CTL
    
    replace_vars $maker_opts genome $gensub_fa
    replace_vars $maker_opts est $RNA
    replace_vars $maker_opts protein $proteins
    # Above three commands as in run_maker1
    replace_vars $maker_opts genome_gff  $gensub.all.gff
    
    replace_vars $maker_opts snaphmm snap.hmm
    replace_vars $maker_opts gmhmm gm_es.mod
    replace_vars $maker_opts augustus_species $augustus_species
    
    replace_vars $maker_opts est2genome 0
    replace_vars $maker_opts protein2genome 0
    
    replace_vars $maker_exe snap \\
    replace_vars $maker_exe augustus \\
    replace_vars 
}

# dummy job to create a file when all previous qsub jobs complete
# see qsmaker.sh for details
finish() {
    touch finished
}

#1. You may wish to copy the evidence files and control files from the last run

#2. Ensure that est2genome and protein2genome is set to 0.

#3. Under the gene predictions section, set augustus_species, snaphmm, gmhmm as needed.

#4. Run MAKER again.

#5. Collect your results as before.


#dir_setup
#sample_fasta
#run_maker1

eval $func

