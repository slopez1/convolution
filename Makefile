.PHONY: serial omp mpi hybrid all clean

serial:
	gcc SERIAL/convolution.c -lm -o SERIAL/serial

mpi:
	mpicc MPI/convolution.c -lm -o MPI/mpi

omp:
	gcc OMP/convolution.c -fopenmp -o OMP/omp

hybrid:
	mpicc HYBRID/convolution.c -lm -fopenmp -o HYBRID/hybrid

all: serial mpi omp hybrid

test_mpi:
	./MPI/run-mpi 01_3 4 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 4 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 4 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 4 im01.ppm  kernel49x49_random.txt

	./MPI/run-mpi 01_3 8 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 8 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 8 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 8 im01.ppm  kernel49x49_random.txt


	./MPI/run-mpi 02_3 4 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 4 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 4 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 4 im02.ppm  kernel49x49_random.txt

	./MPI/run-mpi 02_3 8 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 8 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 8 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 8 im02.ppm  kernel49x49_random.txt


	./MPI/run-mpi 03_3 4 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 4 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 4 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 4 im03.ppm  kernel49x49_random.txt

	./MPI/run-mpi 03_3 8 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 8 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 8 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 8 im03.ppm  kernel49x49_random.txt

	./MPI/run-mpi 04_3 4 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 4 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 4 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 4 im04.ppm  kernel49x49_random.txt

	./MPI/run-mpi 04_3 8 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 8 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 8 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 8 im04.ppm  kernel49x49_random.txt


	./MPI/run-mpi 05_3 4 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 4 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 4 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 4 im05.ppm  kernel49x49_random.txt

	./MPI/run-mpi 05_3 8 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 8 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 8 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 8 im05.ppm  kernel49x49_random.txt

test_serial:
	./SERIAL/run-serial 01_3 im01.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 01_5 im01.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 01_25 im01.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 01_49 im01.ppm  kernel49x49_random.txt

	./SERIAL/run-serial 02_3 im02.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 02_5 im02.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 02_25 im02.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 02_49 im02.ppm  kernel49x49_random.txt

	./SERIAL/run-serial 03_3 im03.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 03_5 im03.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 03_25 im03.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 03_49 im03.ppm  kernel49x49_random.txt

	./SERIAL/run-serial 04_3 im04.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 04_5 im04.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 04_25 im04.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 04_49 im04.ppm  kernel49x49_random.txt

	./SERIAL/run-serial 05_3 im05.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 05_5 im05.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 05_25 im05.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 05_49 im05.ppm  kernel49x49_random.txt

test_omp:
	./OMP/run-omp 01_3 1 im01.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 01_5 1 im01.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 01_25 1 im01.ppm  kernel25x25_random.txt
	./OMP/run-omp 01_49 1 im01.ppm  kernel49x49_random.txt

	./OMP/run-omp 01_3 2 im01.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 01_5 2 im01.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 01_25 2 im01.ppm  kernel25x25_random.txt
	./OMP/run-omp 01_49 2 im01.ppm  kernel49x49_random.txt

	./OMP/run-omp 01_3 4 im01.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 01_5 4 im01.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 01_25 4 im01.ppm  kernel25x25_random.txt
	./OMP/run-omp 01_49 4 im01.ppm  kernel49x49_random.txt

	./OMP/run-omp 02_3 1 im02.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 02_5 1 im02.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 02_25 1 im02.ppm  kernel25x25_random.txt
	./OMP/run-omp 02_49 1 im02.ppm  kernel49x49_random.txt

	./OMP/run-omp 02_3 2 im02.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 02_5 2 im02.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 02_25 2 im02.ppm  kernel25x25_random.txt
	./OMP/run-omp 02_49 2 im02.ppm  kernel49x49_random.txt

	./OMP/run-omp 02_3 4 im02.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 02_5 4 im02.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 02_25 4 im02.ppm  kernel25x25_random.txt
	./OMP/run-omp 02_49 4 im02.ppm  kernel49x49_random.txt

	./OMP/run-omp 03_3 1 im03.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 03_5 1 im03.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 03_25 1 im03.ppm  kernel25x25_random.txt
	./OMP/run-omp 03_49 1 im03.ppm  kernel49x49_random.txt

	./OMP/run-omp 03_3 2 im03.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 03_5 2 im03.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 03_25 2 im03.ppm  kernel25x25_random.txt
	./OMP/run-omp 03_49 2 im03.ppm  kernel49x49_random.txt

	./OMP/run-omp 03_3 4 im03.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 03_5 4 im03.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 03_25 4 im03.ppm  kernel25x25_random.txt
	./OMP/run-omp 03_49 4 im03.ppm  kernel49x49_random.txt

	./OMP/run-omp 04_3 1 im04.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 04_5 1 im04.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 04_25 1 im04.ppm  kernel25x25_random.txt
	./OMP/run-omp 04_49 1 im04.ppm  kernel49x49_random.txt

	./OMP/run-omp 04_3 2 im04.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 04_5 2 im04.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 04_25 2 im04.ppm  kernel25x25_random.txt
	./OMP/run-omp 04_49 2 im04.ppm  kernel49x49_random.txt

	./OMP/run-omp 04_3 4 im04.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 04_5 4 im04.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 04_25 4 im04.ppm  kernel25x25_random.txt
	./OMP/run-omp 04_49 4 im04.ppm  kernel49x49_random.txt

	./OMP/run-omp 05_3 1 im05.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 05_5 1 im05.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 05_25 1 im05.ppm  kernel25x25_random.txt
	./OMP/run-omp 05_49 1 im05.ppm  kernel49x49_random.txt

	./OMP/run-omp 05_3 2 im05.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 05_5 2 im05.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 05_25 2 im05.ppm  kernel25x25_random.txt
	./OMP/run-omp 05_49 2 im05.ppm  kernel49x49_random.txt

	./OMP/run-omp 05_3 4 im05.ppm  kernel3x3_Edge.txt
	./OMP/run-omp 05_5 4 im05.ppm  kernel5x5_Sharpen.txt
	./OMP/run-omp 05_25 4 im05.ppm  kernel25x25_random.txt
	./OMP/run-omp 05_49 4 im05.ppm  kernel49x49_random.txt

test_hybrid:
	./HYBRID/run-hybrid 01_3 1 im01.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 01_5 1 im01.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 01_25 1 im01.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 01_49 1 im01.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 01_3 2 im01.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 01_5 2 im01.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 01_25 2 im01.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 01_49 2 im01.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 01_3 4 im01.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 01_5 4 im01.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 01_25 4 im01.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 01_49 4 im01.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 02_3 1 im02.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 02_5 1 im02.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 02_25 1 im02.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 02_49 1 im02.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 02_3 2 im02.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 02_5 2 im02.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 02_25 2 im02.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 02_49 2 im02.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 02_3 4 im02.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 02_5 4 im02.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 02_25 4 im02.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 02_49 4 im02.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 03_3 1 im03.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 03_5 1 im03.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 03_25 1 im03.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 03_49 1 im03.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 03_3 2 im03.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 03_5 2 im03.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 03_25 2 im03.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 03_49 2 im03.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 03_3 4 im03.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 03_5 4 im03.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 03_25 4 im03.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 03_49 4 im03.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 04_3 1 im04.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 04_5 1 im04.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 04_25 1 im04.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 04_49 1 im04.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 04_3 2 im04.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 04_5 2 im04.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 04_25 2 im04.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 04_49 2 im04.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 04_3 4 im04.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 04_5 4 im04.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 04_25 4 im04.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 04_49 4 im04.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 05_3 1 im05.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 05_5 1 im05.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 05_25 1 im05.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 05_49 1 im05.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 05_3 2 im05.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 05_5 2 im05.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 05_25 2 im05.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 05_49 2 im05.ppm  kernel49x49_random.txt

	./HYBRID/run-hybrid 05_3 4 im05.ppm  kernel3x3_Edge.txt
	./HYBRID/run-hybrid 05_5 4 im05.ppm  kernel5x5_Sharpen.txt
	./HYBRID/run-hybrid 05_25 4 im05.ppm  kernel25x25_random.txt
	./HYBRID/run-hybrid 05_49 4 im05.ppm  kernel49x49_random.txt

clean:
	rm ./SERIAL/serial ./MPI/mpi ./OMP/omp ./HYBRID/hybrid
