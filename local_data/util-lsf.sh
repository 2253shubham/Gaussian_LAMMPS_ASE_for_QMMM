if [ "$HOSTS" ]; then
  ROOT=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
else
  ROOT=$LS_SUBCWD
  WALLTIME=$(bjobs -l $LSB_BATCH_JID | grep -A1 RUNLIMIT | tail -n1 | gawk '{printf "%d", $1*60}')
  NODE=($(uniq $LSB_DJOB_RANKFILE))
  #echo ${NODE[@]} > tmp-$PBS_JOBID
  HOSTS=${NODE[@]}
  nprocpernode=$(echo $LSB_DJOB_NUMPROC ${#NODE[@]} | gawk '{print $1/$2}')
fi
cd "$ROOT"
#######################################################################
FILELOG=log
if [ -z "$nprocperjob" ]; then
    nnodeperjob=$(echo $HOSTS|wc -w) #$PBS_NUM_NODES or $SLURM_NNODES
    nprocperjob=$((nnodeperjob*nprocpernode))
    [ -f ./util.sh ] && source util.sh
fi
HOSTS=${HOSTS// /,}
#######################################################################
echo "Hostname: $(hostname);  ROOT: $ROOT"
echo "Executable: $PROG"
echo "MPIRUN: $MPIRUN"
echo "Hosts: $HOSTS"
echo "nprocpernode: $nprocpernode; nnodeperjob: $nnodeperjob"
echo
