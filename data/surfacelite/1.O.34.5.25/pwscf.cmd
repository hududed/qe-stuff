#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l ncpus=64
#PBS -l mem=128gb
#PBS -l software=espresso
#PBS -l jobfs=1gb
#PBS -l wd
#PBS -q normalbw
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N gc

module load espresso/5.4.0-sandybridge

export TMPDIR=$PBS_JOBFS

mpirun /short/y35/tjf502/QE-5.4.0-bw/bin/pw.x <pwscf.in> pwscf.out
rm -f *.wfc* *.mix*
