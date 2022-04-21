#include <iostream>
#include "solvers.h"
#include "utilities.h"
#include "math.h"

void init_u(double **u, double *x, double *y, int Nx, int Ny) {
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            u[i][j] = signum(-sqrt(x[i]*x[i] + y[j]*y[j]) + 0.1);
        }
    }
}

int main(int argc, char** argv) {
    // Filename
    std::string filename;

    // Parameters for spatial grid
    int Nx, Ny;
    double a, b;
    // Parameters for temporal grid
    double dt, T_max;

    // Receive command line args.
    // If no argument or too little arguments given, the default values 
    // will be set.

    if (argc < 8) {
        a = -1;
        b = 1;
        Nx = 50;
        Ny = 50;
        dt = 0.1;
        T_max = 0.5;
        filename = "data.txt";
    } else {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        a = atof(argv[3]);
        b = atof(argv[4]);
        dt = atof(argv[5]);
        T_max = atof(argv[6]);
        filename = argv[7];

        if (Nx < 2 || Ny < 2) {
            std::cout << "Incorrect 'Nx' or 'Ny' parameters smaller than 2." << std::endl;
            return 1;
        }
        if (dt <= 0) {
            std::cout << "Time increment 'dt' must be positive." << std::endl;
            return 1;
        }
        if (T_max <= 0) {
            std::cout << "Maximum time 'T_max' must be positive." << std::endl;
            return 1;
        }
        if (a >= b) {
            std::cout << "Coordinate 'b' must be larger than 'a'." << std::endl;
            return 1;
        }
        
    }


    // Spatial grid
    //int Nx = 50;
    //int Ny = 50;
    double *x = new double[Nx+1];
    double *y = new double[Ny+1];
    //double a = -1;
    //double b = 1;
    double delta = 0.0001;

    set_grid(x, a, b, Nx);
    set_grid(y, a, b, Ny);

    // Temporal grid
    //double T_max = 0.5;
    //double dt = 0.1;

    // Solution grid - pointer array for column representation of solution u
    double **u = new double*[Nx+1];

    // Allocate u dynamically according to the grid size
    alloc(u, Nx, Ny);

    // Set initial condition for u
    init_u(u, x, y, Nx, Ny);
    
    // Solve the system for the given parameters
    Merson(u, x, y, Nx, Ny, dt, T_max, delta, filename);

    // Delete grid
    delete[] x;
    delete[] y;
    x = nullptr;
    y = nullptr;
    // Delete u - first individual columns, rows after
    dealloc(u, Nx);

    return 0;
}