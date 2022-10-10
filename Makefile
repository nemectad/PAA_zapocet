# GNU c++ compiler
# CC=g++
# MPI compiler
CC=mpic++

# Intel compiler (mainly for profiling)
# IC=icc
# Set Intel environment
#. /opt/intel/oneapi/setvars.sh

SRC=main.cpp solvers.cpp utilities.cpp
FLAGS=-g -std=c++11
#PROFILER=-pg
TRG=main

# Uncomment one
all:
#	$(CC) $(FLAGS) $(PROFILER) $(SRC) -o $(TRG)
#	$(IC) $(SRC) -o $(TRG)
	$(CC) $(FLAGS) $(SRC) -o $(TRG)
#	$(CC) $(SRC) -o $(TRG)

.PHONY : clean
clean:
	-rm main *.o
	-rm -r main.dSYM
