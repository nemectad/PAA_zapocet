#include "solvers.h"
#include "utilities.h"
#include <cmath>
//#define DEBUG


double Laplace_u(double **u, int i, int j, double *x, double *y) {
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];
    return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx) + 
           (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy);
}

double Laplace_u_parallel(double *u, int i, int j, double *x, double *y, int Ny) {
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];

    return (u[(i+1)*Ny + j] - 2*u[i*Ny + j] + u[(i-1)*Ny + j])/(dx*dx) + 
           (u[i*Ny + j+1] - 2*u[i*Ny + j] + u[i*Ny + j-1])/(dy*dy);
}

double F(double **u, int i, int j, double *x, double *y, int Nx, int Ny) {
    // Boundary condition
    if (i == 0 || j == 0 || i == Nx-1 || j == Ny-1) {
        return 0;
    }
    // Laplace
    return Laplace_u(u, i, j, x, y);
}

double F_parallel(double *u, int i, int j, double *x, double *y, int Nx, int Ny) {
    // Boundary condition
    if (i == 0 || j == 0 || i == Nx-1 || j == Ny-1) {
        return 0;
    }
    // Laplace
    return Laplace_u_parallel(u, i, j, x, y, Ny);
}

void update_u(double **u, double ***K, int Nx, int Ny) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = u[i][j] + (K[0][i][j] + 4*K[3][i][j] + K[4][i][j])/6;
        }
    }
}

void update_u_contiguous(double *u, double **K, int arr_len) {
    for (int i = 0; i < arr_len; i++) {
        u[i] = u[i] + (K[0][i] + 4*K[3][i] + K[4][i])/6;
    }
}

void Merson(double **u, double *x, double *y, int Nx, int Ny, double dt, 
            double T, double delta, std::string filename) {
    double t, tau, eps, E, omega;
    bool last = false;
    tau = dt;
    std::ofstream f;
    omega = 0.8;
    t = 0.0;
    
    double **lap_u = new double*[Nx];
    double **lap_k1 = new double*[Nx];
    double **lap_k2 = new double*[Nx];
    double **lap_k3 = new double*[Nx];
    double **lap_k4 = new double*[Nx];
    double **lap_k5 = new double*[Nx];
    double ***K = new double**[5];
    // K matrix initialization - (5)*(Nx)*(Ny) array
    for (int i = 0; i < 5; i++) {
        K[i] = new double*[Nx];
        alloc(K[i], Nx, Ny);
    }
    
    alloc(lap_u, Nx, Ny);
    alloc(lap_k1, Nx, Ny);
    alloc(lap_k2, Nx, Ny);
    alloc(lap_k3, Nx, Ny);
    alloc(lap_k4, Nx, Ny);
    alloc(lap_k5, Nx, Ny);

    // Write initial data:
    f.open(filename, std::ios::out | std::ios_base::app);
    f << Nx << "\t" << Ny << "\n";
    f.close();
    write_data(u, t, x, y, Nx, Ny, filename, &f);

    // Main loop
    while(1) {
        if (std::abs(T-t) < std::abs(tau)) {
            tau = T - t;
            last = true;
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {  
                lap_u[i][j] = F(u, i, j, x, y, Nx, Ny);         
                K[0][i][j] = tau*lap_u[i][j];
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {  
                lap_k1[i][j] = F(K[0], i, j, x, y, Nx, Ny);
                K[1][i][j] = tau*(lap_u[i][j] + lap_k1[i][j]/3);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {  
                lap_k2[i][j] = F(K[1], i, j, x, y, Nx, Ny);
                K[2][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + lap_k2[i][j])/6);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {  
                lap_k3[i][j] = F(K[2], i, j, x, y, Nx, Ny);
                K[3][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + 3*lap_k3[i][j])/8);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {  
                lap_k4[i][j] = F(K[3], i, j, x, y, Nx, Ny);
                K[4][i][j] = tau*(lap_u[i][j] + 0.5*lap_k1[i][j] - 1.5*lap_k3[i][j] +
                             2*lap_k4[i][j]);
            }
        }

        // Compute errors
        E = max_in_mtrx(K, Nx, Ny);
        
        // Update u
        if(E < delta) {
            update_u(u, K, Nx, Ny); 
            t += tau;
            write_data(u, t, x, y, Nx, Ny, filename, &f);
        }
        
        tau = pow(delta/E, 0.2)*omega*tau;

        if (last) 
            break;  
    }
    
    dealloc(lap_k1, Nx);
    dealloc(lap_k2, Nx);
    dealloc(lap_k3, Nx);
    dealloc(lap_k4, Nx);
    dealloc(lap_k5, Nx);
    dealloc(lap_u, Nx);
    for (int i = 0; i < 5; i++) {
        dealloc(K[i], Nx);
    }
    delete[] K;
}

void Merson_parallel(int iproc, MPI_Comm comm, double **u, double *x, double *y, 
                     int Nx, int Ny, double dt, double T, double delta, 
                     std::string filename) {

    double t, tau, eps, E, omega;
    bool last = false;
    tau = dt;
    std::ofstream f;
    omega = 0.8;
    t = 0.0;
    int arr_len = Nx*Ny;


    // Declare local u array -> flattened u
    double u_local[arr_len];

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u_local[i*Ny + j] = u[i][j];
        }
    }
    
    // Allocate buffer matrices
    double *K1 = new double[arr_len];
    double *K2 = new double[arr_len];
    double *K3 = new double[arr_len];
    double *K4 = new double[arr_len];
    double *K5 = new double[arr_len];
    double *lap_u = new double[arr_len];
    double *lap_k1 = new double[arr_len];
    double *lap_k2 = new double[arr_len];
    double *lap_k3 = new double[arr_len];
    double *lap_k4 = new double[arr_len];

    // Root process
    const int root = 0;
    int m1, m2;
    // Separate the main domain into m1*m2 subdomains
    m1 = 2;
    m2 = 2;

    #ifdef DEBUG
    int iter = 0;
    #endif

    // Define the length of a subdomain array
    int stride = (Nx/m1)*(Ny/m2);

    // Declare local array for computation on the subdomain (has the same 
    // number of elements as the corresponding matrix)
    double local_K[stride];
    double local_lap[stride];

    // Write initial data:
    if (iproc == root) {
        f.open(filename, std::ios::out | std::ios_base::app);
        f << Nx << "\t" << Ny << "\n";
        f.close();
        write_contiguous_data(u_local, t, x, y, Nx, Ny, filename, &f);
    }

    // Main loop
    while(1) {
        if (iproc == root) {
            if (std::abs(T-t) < std::abs(tau)) {
                tau = T - t;
                last = true;
            }
        }
        
        // Declare shift integers to move on the computational grid
        int pos_i = (iproc/m1)*(Nx/m1);
        int pos_j = (iproc%m2)*(Ny/m2);
        int pos;

        /**********************
        Compute K1 coefficients
        */

        // Compute the values in each subdomain
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F_parallel(u_local, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                local_K[shift_i+j] = tau*local_lap[shift_i+j];
            }
        }

        /* -----------------
        Gather computed data
        */
        if (iproc == root) {
            gather(local_K, local_lap, lap_u, K1, m1, m2, Nx, Ny, stride, root, comm);
        } else {
            // Other processess do nothing
            MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
        }

        // Broadcasting values to other processes
        MPI_Bcast(K1, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(lap_u, arr_len, MPI_DOUBLE, root, comm);
        
        /**********************
        Compute K2 coefficients
        */
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F_parallel(K1, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                pos = (i+pos_i)*Ny + j + pos_j;
                local_K[shift_i+j] = tau*(lap_u[pos] + local_lap[shift_i+j]/3);
            }
        }

        /* -----------------
        Gather computed data
        */
        if (iproc == root) {
            gather(local_K, local_lap, lap_k1, K2, m1, m2, Nx, Ny, stride, root, comm);
        } else {
            // Other processess do nothing
            MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
        }

        // Broadcasting values to other processes
        MPI_Bcast(K2, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(lap_k1, arr_len, MPI_DOUBLE, root, comm);

        /**********************
        Compute K3 coefficients
        */
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F_parallel(K2, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                pos = (i+pos_i)*Ny + j + pos_j;
                local_K[shift_i+j] = tau*(lap_u[pos] + (lap_k1[pos]+local_lap[shift_i+j])/6);
            }
        }

        /* -----------------
        Gather computed data
        */
        if (iproc == root) {
            gather(local_K, local_lap, lap_k2, K3, m1, m2, Nx, Ny, stride, root, comm);
        } else {
            // Other processess do nothing
            MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
        }

        // Broadcasting values to other processes
        MPI_Bcast(K3, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(lap_k2, arr_len, MPI_DOUBLE, root, comm);


        /**********************
        Compute K4 coefficients
        */
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F_parallel(K3, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                pos = (i+pos_i)*Ny + j + pos_j;
                local_K[shift_i+j] = tau*(lap_u[pos] + (lap_k1[pos]+3*local_lap[shift_i+j])/8);
            }
        }

        /* -----------------
        Gather computed data
        */
        if (iproc == root) {
            gather(local_K, local_lap, lap_k3, K4, m1, m2, Nx, Ny, stride, root, comm);
        } else {
            // Other processess do nothing
            MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
        }

        // Broadcasting values to other processes
        MPI_Bcast(K4, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(lap_k3, arr_len, MPI_DOUBLE, root, comm);

        /**********************
        Compute K5 coefficients
        */
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F_parallel(K4, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                pos = (i+pos_i)*Ny + j + pos_j;
                local_K[shift_i+j] = tau*(lap_u[pos] + 0.5*lap_k1[pos]-1.5*lap_k3[pos] +
                                     2*local_lap[shift_i+j]);
            }
        }

        /* -----------------
        Gather computed data
        */
        if (iproc == root) {
            gather(local_K, local_lap, lap_k4, K5, m1, m2, Nx, Ny, stride, root, comm);
        } else {
            // Other processess do nothing
            MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
        }

        // Broadcasting values to other processes
        MPI_Bcast(K5, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(lap_k4, arr_len, MPI_DOUBLE, root, comm);

        /* ************************
        Compute errors & write data
        */ 
        if (iproc == root) {
            // Declare array of K's
            double *K[5] = {K1, K2, K3, K4, K5};
            E = max_in_arr(K, arr_len);

            // Update solution
            if(E < delta) {
                update_u_contiguous(u_local, K, arr_len); 
                t += tau;
                write_contiguous_data(u_local, t, x, y, Nx, Ny, filename, &f);
            }

            // Compute new value of the timestep
            tau = pow(delta/E, 0.2)*omega*tau;

        } 

        // Broadcast updated values
        MPI_Bcast(u_local, arr_len, MPI_DOUBLE, root, comm);
        MPI_Bcast(&t, 1, MPI_DOUBLE, root, comm);
        MPI_Bcast(&tau, 1, MPI_DOUBLE, root, comm);
        MPI_Bcast(&last, 1, MPI_C_BOOL, root, comm);


        if (last) 
            break;
        
            
        #ifdef DEBUG

        if (iproc == 0) {
            std::cout << "E = " << E << std::endl;
            std::cout << "t = " << t << " tau = " << tau << " last = " << last << std::endl;
        }   
    
        iter++;
        if (iter > 50)
            break;

        #endif
    }

    // Free memory
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] K5;
    delete[] lap_u;
    delete[] lap_k1;
    delete[] lap_k2;
    delete[] lap_k3;
    delete[] lap_k4;

}

void convolution_in_t(int Nx, int Ny, double *x, double *y, double **u, double t) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = convolve(Nx, Ny, x, y, t, x[i], y[j]);
        }
    }
}

double convolve(int Nx, int Ny, double *x, double *y, double t, double m, double n) {
    double conv = 0;
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            conv += Gauss(m-x[i], n-y[j], t)*(signum(-sqrt(x[i]*x[i] + y[j]*y[j]) + 0.1)+1)*dx*dy;
        }
    }
    return conv;
}
