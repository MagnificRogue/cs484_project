SHELL = /bin/bash

CC=g++
MPICC=mpic++

OPT_LEVEL = -O0

BASIC_PERF=$(CC) -g -O0 -std=c++11 -lrt -Wno-write-strings -fpermissive /home/dtubbs2/lib/armadillo-7.800.2/libarmadillo.so.7  

MPI_PERF=$(MPICC) -g -O0 -std=c++11 -lrt -Wno-write-strings -fpermissive /home/dtubbs2/lib/armadillo-7.800.2/libarmadillo.so.7  

all:  main.exe

main.exe: main.cpp
	$(BASIC_PERF) -o main.exe main.cpp GeneticAlgorithm.cpp Constants.h

parallel: main.cpp
	$(BASIC_PERF) -fopenmp -o main_parallel.exe main.cpp GeneticAlgorithm_Parallel.cpp Constants.h


mpi: main_mpi.cpp
	$(MPI_PERF) -fopenmp -o main_parallel_mpi.exe main_mpi.cpp GeneticAlgorithm_Parallel.cpp Constants.h

test: mpi_test.cpp
	$(MPI_PERF) -o test_mpi.exe mpi_test.cpp  Constants.h

clean:
	rm *.exe 
