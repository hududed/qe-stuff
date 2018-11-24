#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l ncpus=56
#PBS -l mem=128gb
#PBS -l software=espresso
#PBS -l jobfs=1gb
#PBS -l wd
#PBS -q normalbw
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N bands

module load espresso/5.4.0-sandybridge

export TMPDIR=$PBS_JOBFS

mpirun /short/y35/tjf502/QE-5.4.0-bw/bin/pw.x <bands.in> bands.out
