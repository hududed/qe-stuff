#PBS -l walltime=48:00:00,nodes=2:ppn=8
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N x

module load OpenMPI/1.8.8-GNU-4.9.3-2.25

cd $PBS_O_WORKDIR


mpirun /home/st/h/hb153250/espresso-5.3.0/bin/pw.x <pwscf.in> pwscf.out

