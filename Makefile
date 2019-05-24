.PHONY: serial omp mpi hybrid all clean

serial:
	gcc SERIAL/convolution.c -lm -o SERIAL/serial

mpi:
	mpicc MPI/convolution.c -lm -o MPI/mpi

all: serial mpi

test_mpi:
	./MPI/run-mpi 01_3 4 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 4 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 4 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 4 im01.ppm  kernel49x49_random.txt
	./MPI/run-mpi 01_99 4 im01.ppm  kernel99x99_random.txt
	./MPI/run-mpi 01_3 8 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 8 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 8 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 8 im01.ppm  kernel49x49_random.txt
	./MPI/run-mpi 01_99 8 im01.ppm  kernel99x99_random.txt
	./MPI/run-mpi 01_3 16 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 16 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 16 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 16 im01.ppm  kernel49x49_random.txt
	./MPI/run-mpi 01_99 16 im01.ppm  kernel99x99_random.txt
	./MPI/run-mpi 01_3 32 im01.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 01_5 32 im01.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 01_25 32 im01.ppm  kernel25x25_random.txt
	./MPI/run-mpi 01_49 32 im01.ppm  kernel49x49_random.txt
	./MPI/run-mpi 01_99 32 im01.ppm  kernel99x99_random.txt


	./MPI/run-mpi 02_3 4 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 4 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 4 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 4 im02.ppm  kernel49x49_random.txt
	./MPI/run-mpi 02_99 4 im02.ppm  kernel99x99_random.txt
	./MPI/run-mpi 02_3 8 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 8 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 8 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 8 im02.ppm  kernel49x49_random.txt
	./MPI/run-mpi 02_99 8 im02.ppm  kernel99x99_random.txt
	./MPI/run-mpi 02_3 16 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 16 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 16 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 16 im02.ppm  kernel49x49_random.txt
	./MPI/run-mpi 02_99 16 im02.ppm  kernel99x99_random.txt
	./MPI/run-mpi 02_3 32 im02.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 02_5 32 im02.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 02_25 32 im02.ppm  kernel25x25_random.txt
	./MPI/run-mpi 02_49 32 im02.ppm  kernel49x49_random.txt
	./MPI/run-mpi 02_99 32 im02.ppm  kernel99x99_random.txt


	./MPI/run-mpi 03_3 4 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 4 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 4 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 4 im03.ppm  kernel49x49_random.txt
	./MPI/run-mpi 03_99 4 im03.ppm  kernel99x99_random.txt
	./MPI/run-mpi 03_3 8 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 8 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 8 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 8 im03.ppm  kernel49x49_random.txt
	./MPI/run-mpi 03_99 8 im03.ppm  kernel99x99_random.txt
	./MPI/run-mpi 03_3 16 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 16 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 16 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 16 im03.ppm  kernel49x49_random.txt
	./MPI/run-mpi 03_99 16 im03.ppm  kernel99x99_random.txt
	./MPI/run-mpi 03_3 32 im03.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 03_5 32 im03.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 03_25 32 im03.ppm  kernel25x25_random.txt
	./MPI/run-mpi 03_49 32 im03.ppm  kernel49x49_random.txt
	./MPI/run-mpi 03_99 32 im03.ppm  kernel99x99_random.txt


	./MPI/run-mpi 04_3 4 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 4 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 4 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 4 im04.ppm  kernel49x49_random.txt
	./MPI/run-mpi 04_99 4 im04.ppm  kernel99x99_random.txt
	./MPI/run-mpi 04_3 8 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 8 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 8 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 8 im04.ppm  kernel49x49_random.txt
	./MPI/run-mpi 04_99 8 im04.ppm  kernel99x99_random.txt
	./MPI/run-mpi 04_3 16 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 16 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 16 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 16 im04.ppm  kernel49x49_random.txt
	./MPI/run-mpi 04_99 16 im04.ppm  kernel99x99_random.txt
	./MPI/run-mpi 04_3 32 im04.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 04_5 32 im04.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 04_25 32 im04.ppm  kernel25x25_random.txt
	./MPI/run-mpi 04_49 32 im04.ppm  kernel49x49_random.txt
	./MPI/run-mpi 04_99 32 im04.ppm  kernel99x99_random.txt


	./MPI/run-mpi 05_3 4 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 4 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 4 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 4 im05.ppm  kernel49x49_random.txt
	./MPI/run-mpi 05_99 4 im05.ppm  kernel99x99_random.txt
	./MPI/run-mpi 05_3 8 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 8 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 8 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 8 im05.ppm  kernel49x49_random.txt
	./MPI/run-mpi 05_99 8 im05.ppm  kernel99x99_random.txt
	./MPI/run-mpi 05_3 16 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 16 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 16 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 16 im05.ppm  kernel49x49_random.txt
	./MPI/run-mpi 05_99 16 im05.ppm  kernel99x99_random.txt
	./MPI/run-mpi 05_3 32 im05.ppm  kernel3x3_Edge.txt
	./MPI/run-mpi 05_5 32 im05.ppm  kernel5x5_Sharpen.txt
	./MPI/run-mpi 05_25 32 im05.ppm  kernel25x25_random.txt
	./MPI/run-mpi 05_49 32 im05.ppm  kernel49x49_random.txt
	./MPI/run-mpi 05_99 32 im05.ppm  kernel99x99_random.txt


#	./MPI/run-mpi 06_3 4 im06.ppm  kernel3x3_Edge.txt
#	./MPI/run-mpi 06_5 4 im06.ppm  kernel5x5_Sharpen.txt
#	./MPI/run-mpi 06_25 4 im06.ppm  kernel25x25_random.txt
#	./MPI/run-mpi 06_49 4 im06.ppm  kernel49x49_random.txt
#	./MPI/run-mpi 06_99 4 im06.ppm  kernel99x99_random.txt


#	./MPI/run-mpi 07_3 4 im07.ppm  kernel3x3_Edge.txt
#	./MPI/run-mpi 07_5 4 im07.ppm  kernel5x5_Sharpen.txt
#	./MPI/run-mpi 07_25 4 im07.ppm  kernel25x25_random.txt
#	./MPI/run-mpi 07_49 4 im07.ppm  kernel49x49_random.txt
#	./MPI/run-mpi 07_99 4 im07.ppm  kernel99x99_random.txt



#	./MPI/run-mpi 10_3 4 im10.ppm  kernel3x3_Edge.txt
#	./MPI/run-mpi 10_5 4 im10.ppm  kkernel5x5_Sharpen.txt
#	./MPI/run-mpi 10_25 4 im10.ppm  kernel25x25_random.txt
#	./MPI/run-mpi 10_49 4 im10.ppm  kernel49x49_random.txt
#	./MPI/run-mpi 10_99 4 im10.ppm  kernel99x99_random.txt


test_serial:
	./SERIAL/run-serial 01_3 im01.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 01_5 im01.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 01_25 im01.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 01_49 im01.ppm  kernel49x49_random.txt
	./SERIAL/run-serial 01_99 im01.ppm  kernel99x99_random.txt


	./SERIAL/run-serial 02_3 im02.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 02_5 im02.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 02_25 im02.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 02_49 im02.ppm  kernel49x49_random.txt
	./SERIAL/run-serial 02_99 im02.ppm  kernel99x99_random.txt


	./SERIAL/run-serial 03_3 im03.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 03_5 im03.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 03_25 im03.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 03_49 im03.ppm  kernel49x49_random.txt
	./SERIAL/run-serial 03_99 im03.ppm  kernel99x99_random.txt


	./SERIAL/run-serial 04_3 im04.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 04_5 im04.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 04_25 im04.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 04_49 im04.ppm  kernel49x49_random.txt
	./SERIAL/run-serial 04_99 im04.ppm  kernel99x99_random.txt


	./SERIAL/run-serial 05_3 im05.ppm  kernel3x3_Edge.txt
	./SERIAL/run-serial 05_5 im05.ppm  kernel5x5_Sharpen.txt
	./SERIAL/run-serial 05_25 im05.ppm  kernel25x25_random.txt
	./SERIAL/run-serial 05_49 im05.ppm  kernel49x49_random.txt
	./SERIAL/run-serial 05_99 im05.ppm  kernel99x99_random.txt


#	./SERIAL/run-serial 06_3 im06.ppm  kernel3x3_Edge.txt
#	./SERIAL/run-serial 06_5 im06.ppm  kernel5x5_Sharpen.txt
#	./SERIAL/run-serial 06_25 im06.ppm  kernel25x25_random.txt
#	./SERIAL/run-serial 06_49 im06.ppm  kernel49x49_random.txt
#	./SERIAL/run-serial 06_99 im06.ppm  kernel99x99_random.txt


#	./SERIAL/run-serial 07_3 im07.ppm  kernel3x3_Edge.txt
#	./SERIAL/run-serial 07_5 im07.ppm  kernel5x5_Sharpen.txt
#	./SERIAL/run-serial 07_25 im07.ppm  kernel25x25_random.txt
#	./SERIAL/run-serial 07_49 im07.ppm  kernel49x49_random.txt
#	./SERIAL/run-serial 07_99 im07.ppm  kernel99x99_random.txt



#	./SERIAL/run-serial 10_3 im10.ppm  kernel3x3_Edge.txt
#	./SERIAL/run-serial 10_5 im10.ppm  kkernel5x5_Sharpen.txt
#	./SERIAL/run-serial 10_25 im10.ppm  kernel25x25_random.txt
#	./SERIAL/run-serial 10_49 im10.ppm  kernel49x49_random.txt
#	./SERIAL/run-serial 10_99 im10.ppm  kernel99x99_random.txt





clean:
	rm serial mpi
