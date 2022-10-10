#include "utilities.h"
#include <cmath>
#include <stdlib.h>

void set_grid(double *arr, double a, double b, int N) {
    double h = (b-a)/(N-1);
    for (int i = 0; i < N; i++) {
        arr[i] = a+h*i;
    }
}

void alloc(double **u, int Ny, int Nx) {
    for (int i = 0; i < Ny; i++) {
        u[i] = new double[Nx];
    }
}

void dealloc(double **u, int Ny) {
    for (int i = 0; i < Ny; i++) {
        delete[] u[i];
    }
    delete[] u;
}

int signum(double x) {
    return ((0 < x) - (x < 0));
}

double Gauss(double x, double y, double t) {
    return exp(-(x*x + y*y)/(4*t))/(4*t*M_PI);
}

double max_in_mtrx(double ***K, int Ny, int Nx) {
    double max = 0;
    double err = 0;
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            err = std::abs(0.2*K[0][i][j] - 0.9*K[2][i][j] + 0.8*K[3][i][j] - 0.1*K[4][i][j])/3;
            max = (max < err) ? err : max;
        }
    }
    return max;
}

double max_in_arr(double **K, int arr_len) {
    double max = 0;
    double err = 0;
    for (int i = 0; i < arr_len; i++) {
        err = std::abs(0.2*K[0][i] - 0.9*K[2][i] + 0.8*K[3][i] - 0.1*K[4][i])/3;
        max = (max < err) ? err : max;
    }
    return max;
}

void write_data(double **u, double t, double *x, double *y, int Ny, 
                int Nx, std::string filename, std::ofstream *f) {
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            *f << y[i] << "\t" << x[j] << "\t" << u[i][j] << "\n";
        }
    }
    f->close();
}

void write_contiguous_data(double *u, double t, double *x, double *y, int Nx, 
                           int Ny, std::string filename, std::ofstream *f) {
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            *f << y[i] << "\t" << x[j] << "\t" << u[i*Ny + j] << "\n";
        }
    }
    f->close();
}

// MODIFY!!!!
void collect_buffers(double *lap, double *K, double *buffer_lap, double *buffer_K, int M1, int M2, int Ny, int Nx) {
    for (int m1 = 0; m1 < M1; m1++) {
        for (int i = 0; i < Nx/M1; i++) {
            for (int m2 = 0; m2 < M2; m2++) {
                for (int j = 0; j < Ny/M2; j++) {
                    lap[i*Ny + j + m2*Ny/M2 + m1*Nx/M1*Ny] = buffer_lap[i*Ny/M2 + j + m1*Ny*(Nx/M1) + m2*(Nx/M1)*(Ny/M2)];
                    K[i*Ny + j + m2*Ny/M2 + m1*Nx/M1*Ny] = buffer_K[i*Ny/M2 + j + m1*Ny*(Nx/M1) + m2*(Nx/M1)*(Ny/M2)];
                }
            }
        }
    }
}

void gather(double *local_K, double *local_lap, double *lap, double *K, int m1, 
            int m2, int count_send, int root, MPI_Comm comm) {
    /**
     * Custom gather function to gather data from all the processes to the root 
     * process.
     * 
     */

    int size = m1*m2;
    int count_recv = count_send;
    // Allocate buffers the size of the computed data
    double *buffer_K = new double[size*count_send];
    double *buffer_lap = new double[size*count_send];

    // Gather all the data into buffers at root
    MPI_Gather(local_K, count_send, MPI_DOUBLE, buffer_K, count_recv, MPI_DOUBLE, root, comm);
    MPI_Gather(local_lap, count_send, MPI_DOUBLE, buffer_lap, count_recv, MPI_DOUBLE, root, comm);
    
    // Collect the flattened buffers
    //collect_buffers(lap, K, buffer_lap, buffer_K, m1, m2, Ny, Nx);

    // Release memory
    delete[] buffer_K;
    delete[] buffer_lap;
}

void share_buffers(MPI_Comm comm, int m1, int m2, int iproc, double *top, 
        double *bottom, double *left, double *right) {

    int M1, M2;
    // Index of the array in the grid
    M1 = iproc/m1;
    M2 = iproc%m2;


}

void set_buffers(double *main_arr, double *left, double *right, double *top, double *bottom, int Ny, int Nx) {
    for (int i = 0; i < Ny; i++) {
        left[i] = main_arr[i*Nx_loc];
        right[i] = main_arr[(i+1)*Nx_loc-1];
    }

    for (int i = 0; i < Nx; i++) {
        top[i] = main_arr[i];
        bottom[i] = main_arr[(Ny_loc-1)*Nx_loc+i];
    }
}