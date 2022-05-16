#include "utilities.h"
#include <cmath>
#include <stdlib.h>

void set_grid(double *arr, double a, double b, int N) {
    double h = (b-a)/(N-1);
    for (int i = 0; i < N; i++) {
        arr[i] = a+h*i;
    }
}

void alloc(double **u, int Nx, int Ny) {
    for (int i = 0; i < Nx; i++) {
        u[i] = new double[Ny];
    }
}

void dealloc(double **u, int Nx) {
    for (int i = 0; i < Nx; i++) {
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

double max_in_mtrx(double ***K, int Nx, int Ny) {
    double max = 0;
    double err = 0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
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

void write_data(double **u, double t, double *x, double *y, int Nx, 
                int Ny, std::string filename, std::ofstream *f) {
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            *f << x[i] << "\t" << y[j] << "\t" << u[i][j] << "\n";
        }
    }
    f->close();
}

void write_contiguous_data(double *u, double t, double *x, double *y, int Nx, 
                           int Ny, std::string filename, std::ofstream *f) {
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            *f << x[i] << "\t" << y[j] << "\t" << u[i*Ny + j] << "\n";
        }
    }
    f->close();
}

/*
void collect_buffers(double *lap, double *K, double *buffer_lap, 
                     double *buffer_K, int m_1, int m_2, int Nx, int Ny) {
    for (int m1 = 0; m1 < m_1; m1++) {
            
        for (int i = 0; i < Nx/m_1; i++) {
            for (int m2 = 0; m2 < m_2; m2++) {
            //int shift_i = i*(Ny/m_2);
                for (int j = 0; j < Ny/m_2; j++) {
                    //lap[(i+shift_k1)*Ny + j+shift_k2] = buffer_lap[shift_i+j];
                    //K[(i+shift_k1)*Ny + j+shift_k2] = buffer_K[shift_i+j];

                    //lap[(i)*Ny + j + m2*Ny/m_2 + m1*(Nx/m_1)*Ny] = buffer_lap[i*(Ny/m_2) + j + m1*(Nx/m_1)*(Ny/m_2) + m2*(Nx)*(Ny/m_2)];
                    //K[(i)*Ny + j + m2*Ny/m_2 + m1*(Nx/m_1)*Ny] = buffer_K[i*(Ny/m_2) + j + m1*(Nx/m_1)*(Ny/m_2) + m2*(Nx)*(Ny/m_2)];

                }
            }
        }
    
    }
}
*/
void collect_buffers(double *lap, double *K, double *buffer_lap, double *buffer_K, int M1, int M2, int Nx, int Ny) {
    for (int m1 = 0; m1 < M1; m1++) {
        for (int i = 0; i < Nx/M1; i++) {
            for (int m2 = 0; m2 < M2; m2++) {
                for (int j = 0; j < Ny/M2; j++) {
                    //std::cout << buffer[i*Ny/M2 + j + m1*Ny*(Nx/M1) + m2*(Nx/M1)*(Ny/M2)] << " ";
                    lap[i*Ny + j + m2*Ny/M2 + m1*Nx/M1*Ny] = buffer_lap[i*Ny/M2 + j + m1*Ny*(Nx/M1) + m2*(Nx/M1)*(Ny/M2)];
                    K[i*Ny + j + m2*Ny/M2 + m1*Nx/M1*Ny] = buffer_K[i*Ny/M2 + j + m1*Ny*(Nx/M1) + m2*(Nx/M1)*(Ny/M2)];
                }
                //std::cout << "\n";
            }
        }
    }
}