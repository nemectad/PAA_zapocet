#ifndef UTILITIES_H
#define UTILITIES_H
#include <fstream>
#include <iostream>

void set_grid(double *arr, double a, double b, int N);
void alloc(double **u, int Nx, int Ny);
void dealloc(double **u, int Nx);
int signum(double x);
double max_in_mtrx(double ***K, int Nx, int Ny);
void write_data(double **u, double t, double *x, double *y, int Nx, 
                int Ny, std::string filename, std::ofstream *f);
double Gauss(double x, double y, double t);

#endif