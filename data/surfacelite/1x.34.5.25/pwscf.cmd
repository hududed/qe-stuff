#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l ncpus=28
#PBS -l mem=128gb
#PBS -l software=espresso
#PBS -l jobfs=1gb
#PBS -q normalbw
#PBS -M hb153250@fh-muenster.de
#PBS -j oe
#PBS -N xdmcooh


export TMPDIR=$PBS_JOBFS


module load openmpi
mpirun /short/y35/tjf502/QE-5.4.0-bw/bin/pw.x <pwscf.in> pwscf.out
rm -f *.wfc* *.mix*
