#ifndef UTILITIES_H
#define UTILITIES_H
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <mpi.h>

void set_grid(double *arr, double a, double b, int N);
void alloc(double **u, int Ny, int Nx);
void dealloc(double **u, int Ny);
int signum(double x);
double max_in_mtrx(double ***K, int Ny, int Nx);
double max_in_arr(double **K, int arr_len);
void write_data(double **u, double t, double *x, double *y, int Ny, 
                int Nx, std::string filename, std::ofstream *f);
void write_contiguous_data(double *u, double t, double *x, double *y, int Nx, 
                           int Ny, std::string filename, std::ofstream *f);
double Gauss(double x, double y, double t);
void collect_buffers(double *lap, double *K, double *buffer_lap, 
                     double *buffer_K, int M1, int M2, int Ny, int Nx);
void gather(double *local_K, double *local_lap, double *lap, double *K, int m1, 
            int m2, int count_send, int root, MPI_Comm comm);
void set_buffers(double *main_arr, double *left, double *right, double *top, 
                 double *bottom, int Ny, int Nx);
void share_buffers(double)

#endif
