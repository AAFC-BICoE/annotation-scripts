#!/bin/bash

# Modifying just the top stanza below in a separate config file should be enough.
#genome_in=/isilon/biodiversity/projects/PRI_Se/assembly/clc-sendo.fasta
#RNA_fasta_in=/isilon/biodiversity/users/cullisj/bug3124/CFIA_contigs.fa
#RNA_gff_in=
#proteins_in=/isilon/biodiversity/users/cullisj/bug3132/B_dendrobatidis/Batde5_best_proteins.fasta
#augustus_species=synchytrium_endobioticum
#cell_type=eukaryote # Bug: no difference if changed. Should change e.g. genemark binary used, among other things.

alt_gm_hmm="../es.mod" # gm model to be used when genemark fails (gm seems to work only with full genome, not subsets).

maker1_dir=maker1
maker_bopts=maker_bopts.ctl
maker_exe=maker_exe.ctl
maker_opts=maker_opts.ctl

augustus1_dir=augustus1
cegma_dir=cegma
snap_dir=snap
genemark_es_dir=genemark_es
genemark_sn_dir=genemark_sn

maker2_dir=maker2
gm_hmm="${genemark_es_dir}/mod/es.mod" # Note: this is a symlink - use readlink to resolve actual path before use.


usage() { echo "Usage: $0 [-w -f <function> | -s <start contig size> -e <end contig size> -f <function>]" 1>&2; exit 1; }
cstart=
cend=
func=
#config_file=
while getopts "wc:s:e:f:" opt; do
    case "${opt}" in
        w)
            contig_range="whole_genome"
            ;;
        c)
            config_file=${OPTARG}
            ;;
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

# Source config file if it exists.
[ -e $config_file ] && source $config_file

config_err() { echo "Error: Config file must define genome_in RNA_fasta_in/RNA_gff_in and proteins_in at a minimum." 1>&2; exit 1; }
echo "Genome in: $genome_in"
echo "RNA fasta in: ${RNA_fasta_in}"
echo "RNA gff in: ${RNA_gff_in}"
echo "Proteins in: $proteins_in"
echo "Augustus species name: $augustus_species"
[[ $genome_in && ($RNA_fasta_in || $RNA_gff_in) && $proteins_in && $augustus_species ]] || config_err

genome_raw=`[ $genome_raw ] || basename $genome_in`
RNA_fasta=` [ ! -z $RNA_fasta_in ] && basename $RNA_fasta_in`
RNA_gff=` [ ! -z $RNA_gff_in ] && basename $RNA_gff_in`
proteins=`[ $proteins ] || basename $proteins_in`

genome="${genome_raw%.*}"
genome_fa=$genome_raw
snap_hmm=$snap_dir/$genome.hmm

# Variables above this line can be redefined in the config.
# Variables below cannot.

if [[ $cstart && $cend ]]; then
    contig_range="$cstart-$cend"
fi

if [[ $contig_range != "whole_genome" ]]; then
    echo "Using contig range $contig_range"
    genome=${genome}_${contig_range}
    genome_fa=$genome.fa
    augustus_species=${augustus_species}_${contig_range} # Used by Maker2
fi

dir_setup()
{
    ln -s $genome_in $genome_raw
    [ ! -z $RNA_fasta_in ] && ln -s $RNA_fasta_in $RNA_fasta
    [ ! -z $RNA_gff_in ] && ln -s $RNA_gff_in $RNA_gff
    ln -s $proteins_in $proteins
}

sample_fasta()
{
    if [[ $contig_range != "whole_genome" ]]; then
        [ -e ./fastaSizes.pl ] || svn export http://biodiversity/svn/source/misc_scripts/fastaSizes.pl
        perl ./fastaSizes.pl -f $genome_raw -r $contig_range -o $genome_fa
    fi
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
    sed -i "s/^$key=[^\s#]*[\s#]/$key=$value #/" $fname
}

run_maker1()
{
    # Run first pass of maker
    mkdir $maker1_dir
    cp $genome_fa $proteins $maker1_dir/
    [ ! -z $RNA_fasta ] && cp $RNA_fasta $maker1_dir/
    [ ! -z $RNA_gff ] && cp $RNA_gff $maker1_dir/
    cd $maker1_dir
    maker -CTL
    replace_vars $maker_opts genome $genome_fa
    [ ! -z $RNA_fasta ] && replace_vars $maker_opts est $RNA_fasta
    [ ! -z $RNA_gff ] && replace_vars $maker_opts est_gff $RNA_gff
    replace_vars $maker_opts protein $proteins
    replace_vars $maker_opts est2genome 1
    #replace_vars $maker_opts protein2genome 1 (only if prokaryotic)
    replace_vars $maker_opts alt_splice 1
    replace_vars $maker_opts TMP \\/state\\/partition1
    replace_vars $maker_bopts blast_type ncbi
    replace_vars $maker_exe RepeatMasker \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/RepeatMasker-open-4-0-3\\/RepeatMasker
    replace_vars $maker_exe exonerate \\/opt\\/bio\\/exonerate\\/bin\\/exonerate
    
    maker_out=$genome.maker.output
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/gff3_merge -d $maker_out/*_master_datastore_index.log
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/fasta_merge -d $maker_out/*_master_datastore_index.log
    sed '/FASTA/q' $genome.all.gff | sed '$d' >$genome.all.nofa.gff
    
    # Line below creates genome.ann, genome.dna files, required by augustus, snap, ..
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker2zff -d $maker_out/*_master_datastore_index.log
    # mv $genome.all* genome.ann genome.dna $maker_out/
    cd ..
}

setup_augustus()
{
    # Only run this if you've never set up before
    mkdir ~/augustus
    cp -r /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/config ~/augustus/
    echo "export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config" >> ~/.bashrc
    echo "export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config" >> ~/.bash_profile
    . ~/.bashrc
}

# Augustus is only for eukaryotes
# Depends on maker1 output
train_augustus1()
{
    mkdir -p augustus1
    cp $maker1_dir/$genome.all.nofa.gff $maker1_dir/genome.dna augustus1/
    cd augustus1
    /isilon/biodiversity/pipelines/maker-2.10/augustus.2.7/scripts/gff2gbSmallDNA.pl $genome.all.nofa.gff genome.dna 1000 genes.gb
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
    cd ..
}

run_cegma()
{
    mkdir -p $cegma_dir
    cd $cegma_dir
    cegma -g ../$genome_fa -o $genome
    # Creates output files $genome.cegma.* *=.dna, .errors, .fa, .gff, .id, .local.gff, .completeness_report 
    # 2. Convert your CEGMA results into SNAP ZFF format by running:
	/isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/cegma2zff ${genome}.cegma.gff ../$genome_fa
    cd ..
}

# Train SNAP using CEGMA output.
train_snap()
{
    mkdir -p $snap_dir
    cp $cegma_dir/genome.ann $cegma_dir/genome.dna $snap_dir/
    cd $snap_dir
    

    # 3. Run all of the following:
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -gene-stats
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -validate
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -categorize 1000
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom uni.ann uni.dna -export 1000 -plus
	mkdir params
	cd params
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/forge ../export.ann ../export.dna
	cd ..
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/hmm-assembler.pl ../$genome params > $genome.hmm
    cd ..
}

# Without running CEGMA first here. This function is deprecated.
# Depends on maker1 output
train_snap_no_cegma()
{
    mkdir -p $snap_dir
    cp $maker1_dir/genome.dna $maker1_dir/genome.ann $snap_dir/
    cd $snap_dir
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -gene-stats
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -validate
    /isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/fathom genome.ann genome.dna -export 1000 -plus
    mkdir params
	cd params
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/forge ../export.ann ../export.dna
	cd ..
	/isilon/biodiversity/pipelines/maker-2.10/snap-2013-16/hmm-assembler.pl ../$genome_fa params > $genome.hmm
	cd ..
}

# Train genemark for eukaryotes
train_genemark_es()
{
    mkdir $genemark_es_dir
    cp $genome_fa $genemark_es_dir/
    #cp $maker1_dir/genome.ann $maker1_dir/genome.dna $genemark_es_dir/
    cd $genemark_es_dir
    
    /isilon/biodiversity/pipelines/maker-2.10/gene-mark-es-2.3e/gmes/gm_es.pl $genome_fa
	# Note: If GeneMark fails, there might be something wrong with your genome.
	# If your contigs are short, try adding --min_contig 10,000 and --max_nnn 5000
    # When completed, the training file is ./mod/es.mod. Note that this is a symlink
    # to another file. If you want to move it, you should just copy the actual file.
    es_mod=$genemark_es_dir/mod/es.mod
    if [ -e $es_mod ]; then
        gm_hmm=`readlink $es_mod`
    elif [ -e $alt_gm_hmm ]; then
        gm_hmm=$alt_gm_hmm
    fi
    cd ..
}

# Train genemark for prokaryotes
train_genemark_sn()
{
    mkdir $genemark_sn_dir
    cp $genome_fa $genemark_sn_dir/
    cd $genemark_sn_dir
    name=gmsn
    gmsn.pl --combine --species $augustus_species -gm --name $name $genome_fa
    # When completed, the training file is <name>_hmm_combined.mod.
    cd ..
}

run_maker2()
{
    mkdir $maker2_dir
    cp $genome_fa $proteins $maker2_dir
    [ ! -z $RNA_fasta ] && cp $RNA_fasta $maker2_dir/
    [ ! -z $RNA_gff ] && cp $RNA_gff $maker2_dir/
    # cp $snap_dir/path/to/snaphmm snap.hmm
    # if [ type=eukaryote ]; then ...
    gm_hmm_rl=`readlink $gm_hmm`
    cp $gm_hmm_rl ${maker2_dir}/gm_es.mod # ln -s to be able to see target here?
    cp $snap_hmm ${maker2_dir}/snap.hmm
    
    cp $maker1_dir/${genome}.all.gff ${maker2_dir}
    cp $maker1_dir/*ctl ${maker2_dir}/
    cd $maker2_dir
    #maker -CTL
    
    replace_vars $maker_opts genome ${genome_fa}
    # replace_vars $maker_opts est $RNA
    [ ! -z $RNA_fasta ] && replace_vars $maker_opts est $RNA_fasta
    [ ! -z $RNA_gff ] && replace_vars $maker_opts est_gff $RNA_gff
    replace_vars $maker_opts protein $proteins
    # Above three commands as in run_maker1
    replace_vars $maker_opts genome_gff  ${genome}.all.gff
    
    replace_vars $maker_opts snaphmm snap.hmm
    replace_vars $maker_opts gmhmm gm_es.mod
    replace_vars $maker_opts augustus_species ${augustus_species}
    
    replace_vars $maker_opts est2genome 0
    replace_vars $maker_opts protein2genome 0
    
    replace_vars $maker_exe snap \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/snap-2013-16\\/snap
    replace_vars $maker_exe augustus \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/augustus.2.7\\/bin\\/augustus
    replace_vars $maker_exe gmhmme3 \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/gene-mark-es-2.3e\\/gmes\\/gmhmme3
    replace_vars $maker_exe probuild \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/gene-mark-es-2.3e\\/gmes\\/probuild
    replace_vars $maker_exe exonerate \\/isilon\\/biodiversity\\/pipelines\\/maker-2.10\\/exonerate-2.2.0\\/bin\\/exonerate
    
    replace_vars $maker_bopts blast_type ncbi
     
    maker_out=${genome}.maker.output
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/gff3_merge -d $maker_out/*_master_datastore_index.log
    /isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/fasta_merge -d $maker_out/*_master_datastore_index.log
    sed '/FASTA/q' $genome.all.gff | sed '$d' >$genome.all.nofa.gff
    
    # Line below creates genome.ann, genome.dna files, required by augustus, snap, ..
    #/isilon/biodiversity/pipelines/maker-2.10/maker-2.10/bin/maker2zff -d $maker_out/*_master_datastore_index.log
    # mv $genome.all* genome.ann genome.dna $maker_out/
    cd ..
}

# dummy job to create a file when all previous qsub jobs complete
# see qsmaker.sh for details
finish() {
    touch finished
}

if [[ $contig_range && $func ]]; then
    eval $func
fi

