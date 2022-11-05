#ifndef SOLVERS_H
#define SOLVERS_H
#include <iostream>
#include <mpi.h>

//double *Laplace_step(double **u, double *x, double *y, int Nx, int Ny);
double Laplace_u(double **u, int i, int j, double *x, double *y);
double Laplace_u_parallel(double *u, int i, int j, double dx, double dy, int Nx);
double F(double **u, int i, int j, double *x, double *y, int Ny, int Nx);
double F_parallel(double *u, double *top, double *bottom, double *left, double *right, 
                  int i, int j, double *x, double *y, int Ny, int Nx, int i_global, int j_global,
                  int Ny_global, int Nx_global);
void Merson(double **u, double *x, double *y, int Ny, int Nx, double dt, 
            double T, double delta, std::string filename);
void Merson_parallel(int iproc, MPI_Comm comm, double **u, double *x, double *y, int Ny, int Nx, double dt, 
            double T, double delta, std::string filename);
void update_u(double **u, double ***K, int Ny, int Nx);
void update_u_contiguous(double *u, double **K, int arr_len);
void convolution_in_t(int Ny, int Nx, double *x, double *y, double **u, double t);
double convolve(int Ny, int Nx, double *x, double *y, double t, double m, double n);

#endif