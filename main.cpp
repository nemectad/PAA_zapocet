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
    // Spatial grid
    int Nx = 50;
    int Ny = 50;
    double *x = new double[Nx+1];
    double *y = new double[Ny+1];
    double a = -1;
    double b = 1;
    double delta = 0.0001;

    set_grid(x, a, b, Nx);
    set_grid(y, a, b, Ny);

    // Temporal grid
    double T_max = 0.5;
    double dt = 0.1;

    // Solution grid - pointer array for column representation of solution u
    double **u = new double*[Nx+1];

    // Allocate u dynamically according to the grid size
    alloc(u, Nx, Ny);

    // Set initial condition for u
    init_u(u, x, y, Nx, Ny);
    
    // Solve the system for the given parameters
    Merson(u, x, y, Nx, Ny, dt, T_max, delta);

    // Delete grid
    delete[] x;
    delete[] y;
    x = nullptr;
    y = nullptr;
    // Delete u - first individual columns, rows after
    dealloc(u, Nx);

    return 0;
}