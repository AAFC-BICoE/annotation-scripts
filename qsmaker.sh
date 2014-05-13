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

# Pass a function name and a set of job ids for the function to hold on.
submit_job() {
    fn=$1
    count=0
    hold_jid_str=
    for j in "$@"; do
        if [ $count ]; then
            hold_jid_str="${hold_jid_str} -hold_jid $j"
        fi
        ((count++))
    done
    # [ $hold_str ] || $hold_str=" -hold_jid 1 "
    qout=`qsub -N "${fn}_${cstart}_${cend}" $hold_jid_str ../qsub_script.sh "../run_maker.sh -s $cstart -e $cend -f $fn"`
    echo $qout 1>&2
    rjid=`echo $qout | awk '{print $3}'`
    echo $rjid
}

# Submit all tasks via qsub.
# Augustus, snap, require output from first maker run
submit_all() {
    dir_setup_jid=`submit_job dir_setup`
    sample_fasta_jid=`submit_job sample_fasta ${dir_setup_jid}`
    genemark_es_jid=`submit_job train_genemark_es ${sample_fasta_jid}`
    maker1_jid=`submit_job run_maker1 ${sample_fasta_jid}`
    augustus1_jid=`submit_job train_augustus1 ${maker1_jid}`
    snap_jid=`submit_job train_snap ${maker1_jid}`
    finished_jid=`submit_job finish ${genemark_es_jid} ${augustus1_jid} ${snap_jid}`
}

if [[ "$func" = "all" ]]; then
    submit_all
else
    qsub -N "${func}_${cstart}_${cend}" ../qsub_script.sh "../run_maker.sh -s $cstart -e $cend -f $func"
fi


