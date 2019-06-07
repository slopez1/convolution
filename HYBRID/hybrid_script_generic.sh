#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Execute the job from the current working directory.
#$ -wd $HOME/convolution/HYBRID/results

## Parallel programming environment (mpich) to instantiate and number of computing slots.
#$ -pe mpich-smp $THREADS

## Passes an enviroment variable to the job
#$ -v OMP_NUM_THREADS=$OMP_THREADS

## The  name  of  the  job.
#$ -N HYBRID_$NAME_$THREADS_$OMP_THREADS

## The folders to save the standard and error outputs.
#$ -o $HOME/convolution/HYBRID/results
#$ -e $HOME/convolution/HYBRID/results

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":"$2}' > $MPICH_MACHINES


## In this line you have to write the command that will execute your application.
mpiexec -f $MPICH_MACHINES -n $NSLOTS $HOME/convolution/HYBRID/hybrid /share/apps/files/convolution/images/$IMG /share/apps/files/convolution/kernel/$KERNEL /dev/null $THREADS

rm -rf $MPICH_MACHINES

