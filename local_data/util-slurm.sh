if [ "$HOSTS" ]; then
  ROOT=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
else
  ROOT=$SLURM_SUBMIT_DIR
  WALLTIME=$(scontrol show job $SLURM_JOB_ID | gawk '/TimeLimit/{seconds=0; split($2,t1,"="); if (split(t1[2],t2,"-")==2) {seconds=t2[1]*24; t2[1]=t2[2];} split(t2[1],t3,":"); print ((seconds+t3[1])*60+t3[2])*60+t3[3]}')
  MEM=$(scontrol show job $SLURM_JOB_ID|gawk '/mem=/{split($1,m1,","); split(m1[2],m2,"="); print m2[2]}'|sed 's/M//')

  # Nodes
  NODE=($(scontrol show hostname "$SLURM_NODELIST"))
  #echo ${NODE[@]} > hosts.txt
  HOSTS=${NODE[@]}

  # CPUs
  hyperthreading=$(scontrol show node ${NODE[0]} | gawk '/ThreadsPerCore=/{split($2,hyperthreading,"="); print hyperthreading[2]}')
  nprocpernode=$((SLURM_CPUS_ON_NODE/hyperthreading))
  ntaskspernode=$SLURM_TASKS_PER_NODE
  [ -z "$SLURM_TASKS_PER_NODE" ] && ntaskspernode=$SLURM_NTASKS_PER_NODE
  ntaskspernode="$(echo $ntaskspernode | sed -E 's/\(x[0-9]+\)//')"
  nproc_requested=$((ntaskspernode*SLURM_CPUS_PER_TASK))
  if [ "$nproc_requested" -lt "$nprocpernode" ]; then
    echo "! == NOTE: Using only $nproc_requested out of $nprocpernode allocated CPUs == !"
    echo
    nprocpernode=$nproc_requested
  fi

  # GPUs
  [ -z "$SLURM_JOB_GPUS" ] && SLURM_JOB_GPUS=$GPU_DEVICE_ORDINAL
  [ -z "$SLURM_JOB_GPUS" ] && SLURM_JOB_GPUS=$SLURM_STEP_GPUS
  IFS=',' read -ra GPUS <<< "$SLURM_JOB_GPUS"
  ngpupernode="${#GPUS[@]}"
fi
cd "$ROOT"
#######################################################################
FILELOG=log
if [ -z "$nprocperjob" ]; then
    nnodeperjob=$(echo $HOSTS|wc -w) #$PBS_NUM_NODES or $SLURM_NNODES
    nprocperjob=$((nnodeperjob*nprocpernode))
    ngpuperjob=$((nnodeperjob*ngpupernode))
    [ -f ./util.sh ] && source util.sh
fi
HOSTS=${HOSTS// /,}
#######################################################################
echo "Hostname: $(hostname);  ROOT: $ROOT"
echo "Executable: $PROG"
echo "MPIRUN: $MPIRUN"
echo "Hosts: $HOSTS"
echo "nnodeperjob: $nnodeperjob; nprocpernode: $nprocpernode; ngpupernode: $ngpupernode"
echo "SLURM_CPUS_ON_NODE: $SLURM_CPUS_ON_NODE; SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"
echo "SLURM_TASKS_PER_NODE: $SLURM_TASKS_PER_NODE; SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo
