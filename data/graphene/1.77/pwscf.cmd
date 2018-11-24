#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l ncpus=16
#PBS -l mem=32gb
#PBS -l software=espresso
#PBS -l jobfs=1gb
#PBS -l wd
#PBS -q normal
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N graphene

module load espresso/5.4.0-sandybridge

export TMPDIR=$PBS_JOBFS

mpirun -np $PBS_NCPUS pw.x <pwscf.in> pwscf.out
rm -f *.wfc* *.mix* *.igk*
