#ifndef SOLVERS_H
#define SOLVERS_H
#include <iostream>

//double *Laplace_step(double **u, double *x, double *y, int Nx, int Ny);
double Laplace_u(double **u, int i, int j, double *x, double *y);
double F(double **u, int i, int j, double *x, double *y, int Nx, int Ny);
void Merson(double **u, double *x, double *y, int Nx, int Ny, double dt, 
            double T, double delta, std::string filename);
void update_u(double **u, double ***K, int Nx, int Ny);

#endif