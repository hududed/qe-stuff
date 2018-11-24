
#!/bin/bash
#PBS -l walltime=15:00:00
#PBS -l ncpus=28
#PBS -l mem=16gb
#PBS -l software=espresso
#PBS -l jobfs=512mb
#PBS -l wd
#PBS -q normalbw
#PBS -M hb153250@fh-muenster.de
#PBS -m ae
#PBS -j oe
#PBS -N xspectra

module load espresso/5.4.0-sandybridge

export TMPDIR=$PBS_JOBFS

mpirun /short/y35/tjf502/QE-5.4.0-bw/bin/xspectra.x <xspectra.in> xspectra.out
rm -f *.wfc* *.mix* *.igk*
