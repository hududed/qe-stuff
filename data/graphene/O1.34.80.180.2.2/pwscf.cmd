#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l ncpus=16
#PBS -l mem=8gb
#PBS -l software=espresso
#PBS -l jobfs=128mb
#PBS -l wd
#PBS -q normal
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N graphene

module load espresso/5.4.0-sandybridge

export TMPDIR=$PBS_JOBFS

mpirun -np $PBS_NCPUS pw.x <pwscf.in> pwscf.out

