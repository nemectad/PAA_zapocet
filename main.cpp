//#include <mpi.h>
#include <iostream>
#include "solvers.h"
#include "utilities.h"
#include <cmath>
#include <stdlib.h>
#include <fstream>

// Uncomment for the analytical model
//#define TEST

void init_u(double **u, double *x, double *y, int Nx, int Ny) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = signum(-sqrt(x[i]*x[i] + y[j]*y[j]) + 0.1)+1;
        }
    }
}

int main(int argc, char* argv[]) {

    // Filename
    std::string filename;

    // Parameters for spatial grid
    int Nx, Ny;
    double a, b;
    // Parameters for temporal grid
    double dt, T_max;
    
    a = -1;
    b = 1;
    Nx = 50;
    Ny = 50;
    dt = 0.02;
    T_max = 0.1;
    filename = "data.txt";

    // Spatial grid
    double *x = new double[Nx];
    double *y = new double[Ny];

    // Error treshold
    double delta = 0.0001;

    set_grid(x, a, b, Nx);
    set_grid(y, a, b, Ny);

    // Solution grid - pointer array for column representation of solution u
    double **u = new double*[Nx];

    // Allocate u dynamically according to the grid size
    alloc(u, Nx, Ny);

    // Set initial condition for u
    init_u(u, x, y, Nx, Ny);
    
    // Solve the system for the given parameters
    Merson_parallel(argc, argv, u, x, y, Nx, Ny, dt, T_max, delta, std::string("parallel_") + filename);


    // Comparison with the analytical model
    #ifdef TEST

    // Computation variables
    int N = 10;
    double *t = new double[N+1];
    init_u(u, x, y, Nx, Ny);

    // Write initial data
    std::ofstream test_f;
    std::string test_filename = std::string("test_") + filename;
    test_f.open(test_filename, std::ios::out | std::ios_base::app);
    test_f << Nx << "\t" << Ny << "\n";
    test_f.close();
    write_data(u, t[0], x, y, Nx, Ny, test_filename, &test_f);

    // Main computation cycle
    for (int i = 1; i <= N; i++) {
        t[i] = t[i-1] + T_max/N;
        convolution_in_t(Nx, Ny, x, y, u, t[i]);
        write_data(u, t[i], x, y, Nx, Ny, test_filename, &test_f);
    }
    delete[] t;
    #endif


    // Delete grid
    delete[] x;
    delete[] y;

    // Delete u - first individual columns, rows after
    dealloc(u, Nx);


    return 0;
}
