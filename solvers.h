#ifndef SOLVERS_H
#define SOLVERS_H
#include <iostream>
#include <mpi.h>

//double *Laplace_step(double **u, double *x, double *y, int Nx, int Ny);
double Laplace_u(double **u, int i, int j, double *x, double *y);
double Laplace_u_parallel(double *u, int i, int j, double *x, double *y, int Ny);
double F(double **u, int i, int j, double *x, double *y, int Nx, int Ny);
double F_parallel(double *u, int i, int j, double *x, double *y, int Nx, int Ny);
void Merson(double **u, double *x, double *y, int Nx, int Ny, double dt, 
            double T, double delta, std::string filename);
void Merson_parallel(int argc, char* argv[], double **u, double *x, double *y, int Nx, int Ny, double dt, 
            double T, double delta, std::string filename);
void update_u(double **u, double ***K, int Nx, int Ny);
void update_u_contiguous(double *u, double **K, int arr_len);
void convolution_in_t(int Nx, int Ny, double *x, double *y, double **u, double t);
double convolve(int Nx, int Ny, double *x, double *y, double t, double m, double n);

#endif