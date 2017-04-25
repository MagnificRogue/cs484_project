SHELL = /bin/bash

CC=g++

OPT_LEVEL = -O0

BASIC_PERF=$(CC) -g -O0 -std=c++11 -lrt -Wno-write-strings -fpermissive /home/dtubbs2/lib/armadillo-7.800.2/libarmadillo.so.7 

all:  main.exe

main.exe: main.cpp
		$(BASIC_PERF) -o main.exe main.cpp GeneticAlgorithm.cpp

parallel: main.cpp
	$(BASIC_PERF) -fopenmp -o main_parallel.exe main.cpp GeneticAlgorithm_Parallel.cpp


clean:
	rm *.exe 
