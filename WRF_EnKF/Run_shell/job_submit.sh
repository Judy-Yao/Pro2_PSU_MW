#!/bin/bash --login
#HOSTTYPE, HOSTPPN are defined in config file
#total_ntasks is defined in run.sh
. $CONFIG_FILE
n=$1  # num of tasks job uses
o=$2  # offset location in total_ntasks, useful for several jobs to run together
ppn=$3  # proc per node for the job
exe=$4  # executable

###stampede
if [[ $HOSTTYPE == "stampede" ]]; then
  export SLURM_TASKS_PER_NODE="$ppn(x$SLURM_NNODES)"
  ibrun -n $n -o $o $exe
  export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
fi

###jet
if [[ $HOSTTYPE == "jet" ]]; then
  if [ $JOB_SUBMIT_MODE == 1 ]; then
  nn=$((($n+$n%$ppn)/$ppn))   #number of nodes used for job = ceiling(n/ppn)
  rm -f nodefile
  for i in `seq 1 $nn`; do 
    cat $PBS_NODEFILE |head -n$((HOSTPPN*($nn*$o/$n+$i))) |tail -n$ppn >> nodefile
  done
  mpiexec.mpirun_rsh -np $n -machinefile nodefile OMP_NUM_THREADS=1 $exe
  #mpiexec -np $n -machinefile nodefile -genv OMP_NUM_THREADS=1 $exe
  #time mpiexec.mpirun_rsh -machinefile \$PBS_NODEFILE -np \$PBS_NP OMP_NUM_THREADS=1
  elif [ $JOB_SUBMIT_MODE == 2 ]; then
    #nodes=`echo "($n+$ppn-1)/$ppn" |bc`
    nodes=`echo "($n+$ppn)/$ppn" |bc`
    jobname=`basename $exe |awk -F. '{print $1}'`
    if [[ $jobname == "wrf" ]] && [[ $wrf_for == "forecast" ]]; then
    nodes=16
    cat << EOF > run_$jobname.sh
#!/bin/bash --login
#SBATCH -A hfip-psu
#SBATCH -J $jobname
#SBATCH -t 8:00:00
#SBATCH -q batch
#SBATCH -p xjet
#SBATCH -N $nodes –ntasks-per-node=$ppn
#SBATCH -ouput=outfile
#SBATCH -error=errfile
cd `pwd`
mpiexec -np $n $exe >& $jobname.log
EOF
    else
    cat << EOF > run_$jobname.sh
#!/bin/bash --login
#SBATCH -A hfip-psu
#SBATCH -J $jobname
#SBATCH -t 2:30:00
#SBATCH -q batch
#SBATCH -p xjet
#SBATCH -N $nodes –ntasks-per-node=$ppn
#SBATCH -ouput=outfile
#SBATCH -error=errfile
cd `pwd`
mpiexec -np $n $exe >& $jobname.log
EOF
    fi
    sbatch run_$jobname.sh >& job_submit.log
    #wait for job to finish
    jobid=`cat job_submit.log |awk -F. '{print $1}'`
    jobstat=1
    until [[ $jobstat == 0 ]]; do
      sleep 1m
      jobstat=`/apps/torque/default/bin/qstat |grep $jobid |awk '{if($5=="R" || $5=="Q") print 1; else print 0;}'`
    done
  fi
fi

###define your own mpiexec here:
