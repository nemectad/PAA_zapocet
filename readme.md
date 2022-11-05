# Diffusion equation solver

**Author**: Tadeáš Němec, 2022

## About
This parallel C++ code uses Merson method to solve 2D diffusion equation using
MPI. 

## How to run
In terminal, move into the code folder. 
Then type the following to compile and run the code:
```shell
make
mpirun -np [n-processes] ./main
```

**Note**: ```n-processes``` has to be an even value.

The file `parallel_data.txt` has been created.

**Note**: To run the program on hardware threads, execute the following command
```shell
mpirun --use-hwthread-cpus -np [n-threads] ./main
```

To visualize the simulation results, run the following Octave/Matlab file:
```shell
octave animate.m
```
The script asks for the file to be plotted, type:
```shell
Enter the name of the data file: parallel_data.txt
```

Then press Enter and the animation starts.

## Benchmarks

Grid size: ```200x200```. *Results without the data writing*.

Serial time: $T_s = 44.520$ s.

Parallel time (4 cores): $T_p = 17.828$ s.

Speedup: $S = T_s/T_p = 2.5$.
