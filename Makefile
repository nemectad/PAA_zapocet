CC=g++
SRC=main.cpp solvers.cpp utilities.cpp
FLAGS=-g
TRG=main

all:
	$(CC) $(FLAGS) $(SRC) -o $(TRG)

.PHONY : clean
clean:
	-rm main *.o
	-rm -r main.dSYM