# Diffusion equation solver

## About
This parallel C++ code uses Merson method to solve 2D diffusion equation using
MPI. 

## How to run
In terminal, move into the code folder. 
Then type the following to compile and run the code:
```shell
make
mpirun -v -np 4 main
```

The file `parallel_data.txt` has been created.
To visualize the simulation results, run the following Octave/Matlab file:
```shell
octave animate.m
```
The script asks for the file to be plotted, type:
```shell
Enter the name of the data file: parallel_data.txt
```

Then press Enter and the animation starts.


Tadeáš Němec, 2022