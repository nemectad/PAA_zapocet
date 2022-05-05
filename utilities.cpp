#include "utilities.h"
#include <cmath>
#include <stdlib.h>

void set_grid(double *arr, double a, double b, int N) {
    double h = (b-a)/N;
    for (int i = 0; i <= N; i++) {
        arr[i] = a+h*i;
    }
}

void alloc(double **u, int Nx, int Ny) {
    for (int i = 0; i <= Nx; i++) {
        u[i] = new double[Ny+1];
    }
}

void dealloc(double **u, int Nx) {
    for (int i = 0; i <= Nx; i++) {
        delete[] u[i];
    }
    delete[] u;
    //u = nullptr;
}

int signum(double x) {
    return ((0 < x) - (x < 0));
}

double Gauss(double x, double y, double t) {
    return exp(-(x*x + y*y)/(4*t))/(4*t*M_PI);
}

double max_in_mtrx(double ***K, int Nx, int Ny) {
    double max = 0;
    double err = 0;
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            err = std::abs(0.2*K[0][i][j] - 0.9*K[2][i][j] + 0.8*K[3][i][j] - 0.1*K[4][i][j])/3;
            max = (max < err) ? err : max;
        }
    }
    return max;
}

void write_data(double **u, double t, double *x, double *y, int Nx, 
                int Ny, std::string filename, std::ofstream *f) {
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            *f << x[i] << "\t" << y[j] << "\t" << u[i][j] << "\n";
        }
    }
    f->close();
}
