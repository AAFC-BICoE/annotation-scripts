usage() { echo "Usage: $0 [-w -f <function> | -s <start contig size> -e <end contig size> -f <function>" 1>&2; exit 1; }

qsub_script=$HOME/svn/Maker/qsub_script.sh
cstart=
cend=
func=
config_args=
do_submit_all=
do_submit_all_from=
while getopts "awc:s:e:f:r:" opt; do
    case "${opt}" in
        a)
            do_submit_all=1
            ;;
        r)
            do_submit_all_from=${OPTARG}
            ;;
        w)
            contig_range="whole_genome"
            ;;
        c)
            config_args=" -c ${OPTARG} "
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

if [[ $cstart && $cend && $func ]]; then
    contig_range="$cstart-$cend"
    contig_args=" -s $cstart -e $cend "
elif [[ $contig_range && $func ]]; then
    contig_args=" -w "
else
    #usage
    x=
fi

# Pass a function name and a set of job ids for the function to hold on.
submit_job() {
    fn=$1
    count=0
    hold_jid_str=
    for j in "$@"; do
        if [ $count -gt 0 ]; then
            hold_jid_str="${hold_jid_str} -hold_jid $j"
        fi
        ((count++))
    done
    # [ $hold_str ] || $hold_str=" -hold_jid 1 "
    qscmd="qsub -N ${fn}_${contig_range} $hold_jid_str $qsub_script \"run_maker.sh $contig_args $config_args -f $fn\""
    echo $qscmd 1>&2
    #qout=`$qscmd`
    qout=`qsub -N ${fn}_${contig_range} $hold_jid_str $qsub_script "run_maker.sh $contig_args $config_args -f $fn"`
    echo $qout 1>&2
    #rjid=`echo $qout | awk '{print $3}'`
    rjid=`echo $qout | perl -ne 'if (/job\s([0-9]+)\s/) { print $1; }'`
    echo $rjid
}

# Submit all tasks via qsub.
# Augustus, snap, require output from first maker run
submit_all() {
    dir_setup_jid=`submit_job dir_setup`
    sample_fasta_jid=`submit_job sample_fasta ${dir_setup_jid}`
    genemark_es_jid=`submit_job train_genemark_es ${sample_fasta_jid}`
    cegma_jid=`submit_job run_cegma ${sample_fasta_jid}`
    maker1_jid=`submit_job run_maker1 ${sample_fasta_jid}`
    augustus1_jid=`submit_job train_augustus1 ${maker1_jid}`
    snap_jid=`submit_job train_snap ${maker1_jid} ${cegma_jid}`
    maker2_jid=`submit_job run_maker2 ${augustus1_jid} ${snap_jid} ${genemark_es_jid}`
    finished_jid=`submit_job finish ${maker2_jid}`
}

# Submit all tasks beyond a given step number via qsub.
submit_all_from() {
    dir_setup_jid=1; sample_fasta_jid=1; genemark_es_jid=1; cegma_jid=1; maker1_jid=1;
    augustus1_jid=1; snap_jid=1; maker2_jid=1; finished_jid=1;
    case "$1" in 
        'dir_setup')
            dir_setup_jid=`submit_job dir_setup`
            ;&
        'sample_fasta')
            sample_fasta_jid=`submit_job sample_fasta ${dir_setup_jid}`
            ;&
        'genemark_es')
            genemark_es_jid=`submit_job train_genemark_es ${sample_fasta_jid}`
            ;&
        'cegma')
            cegma_jid=`submit_job run_cegma ${sample_fasta_jid}`
            ;&
        'maker1')
            maker1_jid=`submit_job run_maker1 ${sample_fasta_jid}`
            ;&
        'augustus1')
            augustus1_jid=`submit_job train_augustus1 ${maker1_jid}`
            ;&
        'snap')
            snap_jid=`submit_job train_snap ${maker1_jid} ${cegma_jid}`
            ;&
        'maker2')
            maker2_jid=`submit_job run_maker2 ${augustus1_jid} ${snap_jid} ${genemark_es_jid}`
            ;&
        'finished')
            finished_jid=`submit_job finish ${maker2_jid}`
            ;;
    esac
}

if [[ ! -z $do_submit_all ]]; then
    submit_all
elif [[ ! -z $do_submit_all_from ]]; then
    submit_all_from $do_submit_all_from
else
    qscmd="qsub -N ${func}_${contig_range} $qsub_script \"run_maker.sh $contig_args $config_args -f $func\""
    echo $qscmd 1>&2
    eval $qscmd
    #qsub -N ${func}_${contig_range} $qsub_script "run_maker.sh $contig_args $config_args -f $func"
fi


