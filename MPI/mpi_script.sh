#!/bin/bash

## Specifies the interpreting shell for the job.
#$ -S /bin/bash

## Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -V

## Execute the job from the current working directory.
#$ -wd $HOME/convolution/MPI/results

## Parallel programming environment (mpich) to instantiate and number of computing slots.
#$ -pe mpich 8

## The  name  of  the  job.
#$ -N MPI_05_49_8

## The folders to save the standard and error outputs.
#$ -o $HOME/convolution/MPI/results
#$ -e $HOME/convolution/MPI/results

MPICH_MACHINES=$TMPDIR/mpich_machines
cat $PE_HOSTFILE | awk '{print $1":"$2}' > $MPICH_MACHINES


## In this line you have to write the command that will execute your application.
mpiexec -f $MPICH_MACHINES -n $NSLOTS $HOME/convolution/MPI/mpi /share/apps/files/convolution/images/im05.ppm /share/apps/files/convolution/kernel/kernel49x49_random.txt /dev/null 8

rm -rf $MPICH_MACHINES

