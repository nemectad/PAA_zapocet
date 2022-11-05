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

void share_buffers(MPI_Comm comm, int m1, int m2, int M1, int M2, int iproc, double *top, 
        double *bottom, double *left, double *right, int Ny_loc, int Nx_loc) {

    int proc_receiver;
    int proc_sender;

    // set temporary buffers
    double right_temp[Ny_loc];
    double top_temp[Nx_loc];

    for (int i = 0; i < Ny_loc; i++) {
        right_temp[i] = right[i];
    }
    for (int i = 0; i < Nx_loc; i++) {
        top_temp[i] = top[i];
    }

    
    // iproc bottom to top of below    
    if (M1 + 1 < m1) {
        proc_receiver = iproc + m2;
        MPI_Send(bottom, Nx_loc, MPI_DOUBLE, proc_receiver, 0, comm);
    }
    if (M1 != 0) {
        proc_sender = iproc - m2;
        MPI_Recv(top, Nx_loc, MPI_DOUBLE, proc_sender, 0, comm, MPI_STATUS_IGNORE);
    }
    
    // iproc top to bottom of up
    if (M1 - 1 >= 0) {
        proc_receiver = iproc - m2;
        MPI_Send(top_temp, Nx_loc, MPI_DOUBLE, proc_receiver, 0, comm);
    }
    if (M1 != m1 - 1) {
        proc_sender = iproc + m2;
        MPI_Recv(bottom, Nx_loc, MPI_DOUBLE, proc_sender, 0, comm, MPI_STATUS_IGNORE);
    }
    
    // iproc left to right of left
    if (M2 - 1 >= 0) {
        proc_receiver = iproc - 1;
        MPI_Send(left, Ny_loc, MPI_DOUBLE, proc_receiver, 0, comm);
    }
    if (M2 != m2 - 1) {
        proc_sender = iproc + 1;
        MPI_Recv(right, Ny_loc, MPI_DOUBLE, proc_sender, 0, comm, MPI_STATUS_IGNORE);
    }

    // iproc right to left of right
    if (M2 + 1 < m2) {
        proc_receiver = iproc + 1;
        MPI_Send(right_temp, Ny_loc, MPI_DOUBLE, proc_receiver, 0, comm);
    }
    if (M2 != 0) {
        proc_sender = iproc - 1;
        MPI_Recv(left, Ny_loc, MPI_DOUBLE, proc_sender, 0, comm, MPI_STATUS_IGNORE);
    }
    
    //MPI_Barrier(comm);
}

void set_buffers(double *main_arr, double *left, double *right, double *top, double *bottom, int Ny, int Nx) {
    for (int i = 0; i < Ny; i++) {
        left[i] = main_arr[i*Nx];
        right[i] = main_arr[(i+1)*Nx-1];
    }

    for (int i = 0; i < Nx; i++) {
        top[i] = main_arr[i];
        bottom[i] = main_arr[(Ny-1)*Nx+i];
    }
}

void collect_and_write_u(MPI_Comm comm, double *local_u, int m1, int m2, 
                         int count_send, double t, double *x, double *y, int Ny, 
                         int Nx, int root, std::string filename, std::ofstream *f) {
    int size = m1*m2;
    int Nx_loc = Nx/m2;
    int Ny_loc = Ny/m1;
    int count_recv = count_send;
    // Allocate buffers the size of the computed data
    double *u = new double[size*count_send];

    // Gather all the data into buffers at root
    MPI_Gather(local_u, count_send, MPI_DOUBLE, u, count_recv, MPI_DOUBLE, root, comm);
    
    f->open(filename, std::ios::out | std::ios_base::app);
    *f << t << "\n";

    for (int M1 = 0; M1 < m1; M1++) {
        for (int i = 0; i < Ny_loc; i++) {
            for (int M2 = 0; M2 < m2; M2++) {
                for (int j = 0; j < Nx_loc; j++) {
                    *f << y[i+M1*Ny_loc] << "\t" << x[j+M2*Nx_loc] << "\t" << 
                        u[M1*Ny_loc*Nx + i*Nx_loc + M2*Nx_loc*Ny_loc + j] << "\n";
                }
            }
        }
    }

    f->close();

    delete [] u;
}
