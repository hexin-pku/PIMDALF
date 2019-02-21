#PBS -N mqcivr
#PBS -l nodes=5:ppn=12,walltime=15:00:00
#PBS -S /bin/bash
#PBS -q default
echo $PBS_NODEFILE
echo `cat $PBS_NODEFILE`

export GMPICONF=nodeinfo/PBS_JOBID

cd $PBS_O_WORKDIR

time mpirun -np 60 ./dyn.x >dyn.out

