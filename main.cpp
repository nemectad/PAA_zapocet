#include <iostream>
#include "solvers.h"
#include "utilities.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <ctime>
#include <mpi.h>

// Uncomment for the analytical model
//#define TEST
// Uncomment for serial code
//#define SERIAL
// Uncomment for parallel code
#define PARALLEL

void init_u(double **u, double *x, double *y, int Ny, int Nx) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            u[i][j] = signum(-sqrt(y[i]*y[i] + x[j]*x[j]) + 0.1)+1;
        }
    }
}

int main(int argc, char* argv[]) {

    // Filename
    std::string filename;

    // Parameters for spatial grid
    int Nx, Ny;
    double a, b, dx, dy;
    // Parameters for temporal grid
    double dt, T_max;
    
    a = -1;
    b = 1;
    Nx = 50;
    Ny = 50;
    dt = 0.02;
    T_max = 0.3;
    filename = "data.txt";

    // Spatial grid
    double *x = new double[Nx];
    double *y = new double[Ny];

    // Solution grid - pointer array for column representation of solution u
    double **u = new double*[Ny];

    // Error treshold
    double delta = 0.0001;

    set_grid(x, a, b, Nx);
    set_grid(y, a, b, Ny);

    // Allocate u dynamically according to the grid size
    alloc(u, Ny, Nx);

    

    // Set initial condition for u
    init_u(u, x, y, Ny, Nx);
    
    #ifdef PARALLEL
    // ****************************************************
    // Solve the system for the given parameters - parallel

    int nproc, iproc;
        
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    if(nproc%2 != 0) {
        std::cout << "Program is meant to run with odd number of processes!" << std::endl;
    }
    
    Merson_parallel(iproc, MPI_COMM_WORLD, u, x, y, Ny, Nx, dt, T_max, delta, std::string("parallel_") + filename);

    MPI_Finalize();
    #endif

  
    #ifdef SERIAL
    // ****************************************************
    // Solve the system for the given parameters - serial

    // Set initial condition for u
    //init_u(u, x, y, Ny, Nx);

    // **************************************************
    // Solve the system for the given parameters - serial
    
    Merson(u, x, y, Ny, Nx, dt, T_max, delta, std::string("serial_") + filename);

    #endif

    // Comparison with the analytical model
    #ifdef TEST

    // Computation variables
    int N = 10;
    double t = 0;
    dt = T_max/N;
    init_u(u, x, y, Nx, Ny);

    // Write initial data
    std::ofstream test_f;
    std::string test_filename = std::string("test_") + filename;
    test_f.open(test_filename, std::ios::out | std::ios_base::app);
    test_f << Nx << "\t" << Ny << "\n";
    test_f.close();
    write_data(u, t, x, y, Nx, Ny, test_filename, &test_f);

    // Main computation cycle
    for (int i = 1; i <= N; i++) {
        t = t + dt;
        convolution_in_t(Nx, Ny, x, y, u, t);
        write_data(u, t, x, y, Nx, Ny, test_filename, &test_f);
    }
    
    #endif


    // Delete grid
    delete[] x;
    delete[] y;

    // Delete u - first individual columns, rows after
    dealloc(u, Ny);


    return 0;
}
