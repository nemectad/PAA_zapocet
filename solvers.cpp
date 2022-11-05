#include "solvers.h"
#include "utilities.h"
#include <cmath>

//#define DEBUG


double Laplace_u(double **u, int i, int j, double *x, double *y) {
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];
    return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dy*dy) + 
           (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dx*dx);
}

double Laplace_u_parallel(double *u, int i, int j, double dx, double dy, int Nx) {
    return (u[(i+1)*Nx + j] - 2*u[i*Nx + j] + u[(i-1)*Nx + j])/(dx*dx) + 
           (u[i*Nx + j+1] - 2*u[i*Nx + j] + u[i*Nx + j-1])/(dy*dy);
}

double F(double **u, int i, int j, double *x, double *y, int Ny, int Nx) {
    // Boundary condition
    if (i == 0 || j == 0 || i == Ny-1 || j == Nx-1) {
        return 0;
    }
    // Laplace
    return Laplace_u(u, i, j, x, y);
}

double F_parallel(double *u, double *top, double *bottom, double *left, double *right,
                  int i, int j, double dx, double dy, int Ny, int Nx, int i_global, int j_global,
                  int Ny_global, int Nx_global) {
    double res = 0;
    bool x, y;
    x = false;
    y = false;

    if (i_global == 0 || i_global == Ny_global-1) return 0.;
    if (j_global == 0 || j_global == Nx_global-1) return 0.;
    
    if (j == 0) {
        // Left boundary
        res += (u[i*Nx + j+1] - 2*u[i*Nx + j] + left[i])/(dx*dx);
        x = true;
    } else if (j == Nx-1) {
        // Right boundary
        res += (right[i] - 2*u[i*Nx + j] + u[i*Nx + j-1])/(dx*dx);
        x = true;
    }
    
    if (i == 0) {
        // Top boundary
        res += (u[(i+1)*Nx + j] - 2*u[i*Nx + j] + top[j])/(dy*dy);
        y = true;
    } else if (i == Ny-1) {
        // Bottom boundary
        res += (bottom[j] - 2*u[i*Nx + j] + u[(i-1)*Nx + j])/(dy*dy);
        y = true;
    }

    if (!y) {
        res += (u[(i+1)*Nx + j] - 2*u[i*Nx + j] + u[(i-1)*Nx + j])/(dy*dy);
    } 
    if (!x) {
        res += (u[i*Nx + j+1] - 2*u[i*Nx + j] + u[i*Nx + j-1])/(dx*dx);
    }
    return res;
    
}

void update_u(double **u, double ***K, int Ny, int Nx) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            u[i][j] = u[i][j] + (K[0][i][j] + 4*K[3][i][j] + K[4][i][j])/6;
        }
    }
}

void update_u_contiguous(double *u, double **K, int arr_len) {
    for (int i = 0; i < arr_len; i++) {
        u[i] = u[i] + (K[0][i] + 4*K[3][i] + K[4][i])/6;
    }
}

void Merson(double **u, double *x, double *y, int Ny, int Nx, double dt, 
            double T, double delta, std::string filename) {
    double t, tau, eps, E, omega;
    bool last = false;
    tau = dt;
    std::ofstream f;
    omega = 0.8;
    t = 0.0;
    
    double **lap_u = new double*[Ny];
    double **lap_k1 = new double*[Ny];
    double **lap_k2 = new double*[Ny];
    double **lap_k3 = new double*[Ny];
    double **lap_k4 = new double*[Ny];
    double **lap_k5 = new double*[Ny];
    double ***K = new double**[5];
    // K matrix initialization - (5)*(Nx)*(Ny) array
    for (int i = 0; i < 5; i++) {
        K[i] = new double*[Ny];
        alloc(K[i], Ny, Nx);
    }
    
    alloc(lap_u, Ny, Nx);
    alloc(lap_k1, Ny, Nx);
    alloc(lap_k2, Ny, Nx);
    alloc(lap_k3, Ny, Nx);
    alloc(lap_k4, Ny, Nx);
    alloc(lap_k5, Ny, Nx);

    // Write initial data:
    f.open(filename, std::ios::out | std::ios_base::app);
    f << Ny << "\t" << Nx << "\n";
    f.close();
    write_data(u, t, x, y, Ny, Nx, filename, &f);

    // Main loop
    while(1) {
        if (std::abs(T-t) < std::abs(tau)) {
            tau = T - t;
            last = true;
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {  
                lap_u[i][j] = F(u, i, j, x, y, Ny, Nx);         
                K[0][i][j] = tau*lap_u[i][j];
            }
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {  
                lap_k1[i][j] = F(K[0], i, j, x, y, Ny, Nx);
                K[1][i][j] = tau*(lap_u[i][j] + lap_k1[i][j]/3);
            }
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {  
                lap_k2[i][j] = F(K[1], i, j, x, y, Ny, Nx);
                K[2][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + lap_k2[i][j])/6);
            }
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {  
                lap_k3[i][j] = F(K[2], i, j, x, y, Ny, Nx);
                K[3][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + 3*lap_k3[i][j])/8);
            }
        }

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {  
                lap_k4[i][j] = F(K[3], i, j, x, y, Ny, Nx);
                K[4][i][j] = tau*(lap_u[i][j] + 0.5*lap_k1[i][j] - 1.5*lap_k3[i][j] +
                             2*lap_k4[i][j]);
            }
        }

        // Compute errors
        E = max_in_mtrx(K, Ny, Nx);
        
        // Update u
        if(E < delta) {
            update_u(u, K, Ny, Nx); 
            t += tau;
            write_data(u, t, x, y, Ny, Nx, filename, &f);
        }
        
        tau = pow(delta/E, 0.2)*omega*tau;

        if (last) 
            break;  
    }
    
    dealloc(lap_k1, Ny);
    dealloc(lap_k2, Ny);
    dealloc(lap_k3, Ny);
    dealloc(lap_k4, Ny);
    dealloc(lap_k5, Ny);
    dealloc(lap_u, Ny);
    for (int i = 0; i < 5; i++) {
        dealloc(K[i], Ny);
    }
    delete[] K;

    std::cout << "DONE" << std::endl;
}

void Merson_parallel(int iproc, MPI_Comm comm, double **u, double *x, double *y, 
                     int Ny, int Nx, double dt, double T, double delta, 
                     std::string filename) {

    double t, tau, eps, E, E_global, omega, dx, dy;
    bool last = false;
    tau = dt;
    std::ofstream f;
    omega = 0.8;
    t = 0.0;
    int arr_len = Ny*Nx;
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Root process
    const int root = 0;
    int m1, m2;
    int Ny_loc, Nx_loc;
    
    // Separate the main domain into m1*m2 subdomains
    m1 = 2; // # of rows
    m2 = size/2; // # of columns
    Ny_loc = Ny/m1;
    Nx_loc = Nx/m2;

    // Local array length
    int loc_arr_len = Ny_loc*Nx_loc;

    dx = x[1]-x[0];
    dy = y[1]-y[0];

    // Allocate buffer matrices
    double *K1 = new double[loc_arr_len];
    double *K2 = new double[loc_arr_len];
    double *K3 = new double[loc_arr_len];
    double *K4 = new double[loc_arr_len];
    double *K5 = new double[loc_arr_len];
    double *lap_u = new double[loc_arr_len];
    double *lap_k1 = new double[loc_arr_len];
    double *lap_k2 = new double[loc_arr_len];
    double *lap_k3 = new double[loc_arr_len];
    double *lap_k4 = new double[loc_arr_len];

    // Declare shift integers to move on the computational grid
    int pos_i = (iproc/m2)*Ny_loc;
    int pos_j = (iproc%m2)*Nx_loc;
    int pos;
    int shift_i;

    // Declare local u array -> flattened u
    double local_u[loc_arr_len];

    // Declare local array for computation on the subdomain (has the same 
    // number of elements as the corresponding matrix)
    double local_K[loc_arr_len];
    double local_lap[loc_arr_len];
    double *top = new double[Nx_loc];
    double *bottom = new double[Nx_loc];
    double *left = new double[Ny_loc];
    double *right = new double[Ny_loc];

    int M1, M2;
    // Index of the array in the grid
    M1 = iproc/m2;
    M2 = iproc%m2;


    for (int i = 0; i < Ny_loc; i++) {
        for (int j = 0; j < Nx_loc; j++) {
            local_u[i*Nx_loc + j] = u[i+pos_i][j+pos_j];
        }
    }

    #ifdef DEBUG
    int iter = 0;
    #endif


    // Write initial data:
    if (iproc == root) {
        f.open(filename, std::ios::out | std::ios_base::app);
        f << Ny << "\t" << Nx << "\n";
        f.close();
        write_data(u, t, x, y, Ny, Nx, filename, &f);
    }

    // Main loop
    while(1) {
        if (std::abs(T-t) < std::abs(tau)) {
            tau = T - t;
            last = true;
        }
        

        /**********************
         * Set top, bottom, right and left buffers
        */
        
        set_buffers(local_u, left, right, top, bottom, Ny_loc, Nx_loc);

        /**********************
         * Share buffers 
        */
        share_buffers(comm, m1, m2, M1, M2, iproc, top, bottom, left, right, Ny_loc, Nx_loc);

        /**********************
        Compute K1 coefficients
        */

        // Compute the values in each subdomain & set buffers
        for (int i = 0; i < Ny_loc; i++) {
            shift_i = i*Nx_loc;
            for (int j = 0; j < Nx_loc; j++) {  
                lap_u[shift_i+j] = F_parallel(local_u, top, bottom, left, right, i, j, dx, dy, Ny_loc, Nx_loc, i + pos_i, j + pos_j, Ny, Nx);         
                K1[shift_i+j] = tau*lap_u[shift_i+j];
            }
        }

        /**********************
         * Set top, bottom, right and left buffers
        */
        set_buffers(K1, left, right, top, bottom, Ny_loc, Nx_loc);

        /**********************
         * Share buffers 
        */
        share_buffers(comm, m1, m2, M1, M2, iproc, top, bottom, left, right, Ny_loc, Nx_loc);

        /**********************
        Compute K2 coefficients
        */
        for (int i = 0; i < Ny_loc; i++) {
            shift_i = i*Nx_loc;
            for (int j = 0; j < Nx_loc; j++) {  
                lap_k1[shift_i+j] = F_parallel(K1, top, bottom, left, right, i, j, dx, dy, Ny_loc, Nx_loc, i + pos_i, j + pos_j, Ny, Nx);
                K2[shift_i+j] = tau*(lap_u[shift_i+j] + lap_k1[shift_i+j]/3);
            }
        }
        
        /**********************
         * Set top, bottom, right and left buffers
        */
        set_buffers(K2, left, right, top, bottom, Ny_loc, Nx_loc);

        /**********************
         * Share buffers 
        */
        share_buffers(comm, m1, m2, M1, M2, iproc, top, bottom, left, right, Ny_loc, Nx_loc);

        /**********************
        Compute K3 coefficients
        */
        for (int i = 0; i < Ny_loc; i++) {
            shift_i = i*Nx_loc;
            for (int j = 0; j < Nx_loc; j++) {  
                lap_k2[shift_i+j] = F_parallel(K2, top, bottom, left, right, i, j, dx, dy, Ny_loc, Nx_loc, i + pos_i, j + pos_j, Ny, Nx);
                K3[shift_i+j] = tau*(lap_u[shift_i+j] + (lap_k1[shift_i+j]+lap_k2[shift_i+j])/6);
            }
        }

        /**********************
         * Set top, bottom, right and left buffers
        */
        set_buffers(K3, left, right, top, bottom, Ny_loc, Nx_loc);

        /**********************
         * Share buffers 
        */
        share_buffers(comm, m1, m2, M1, M2, iproc, top, bottom, left, right, Ny_loc, Nx_loc);

        /**********************
        Compute K4 coefficients
        */
        for (int i = 0; i < Ny_loc; i++) {
            shift_i = i*Nx_loc;
            for (int j = 0; j < Nx_loc; j++) {  
                lap_k3[shift_i+j] = F_parallel(K3, top, bottom, left, right, i, j, dx, dy, Ny_loc, Nx_loc, i + pos_i, j + pos_j, Ny, Nx);
                K4[shift_i+j] = tau*(lap_u[shift_i+j] + (lap_k1[shift_i+j]+3*lap_k3[shift_i+j])/8);
            }
        }

        /**********************
         * Set top, bottom, right and left buffers
        */
        set_buffers(K4, left, right, top, bottom, Ny_loc, Nx_loc);

        /**********************
         * Share buffers 
        */
        share_buffers(comm, m1, m2, M1, M2, iproc, top, bottom, left, right, Ny_loc, Nx_loc);

        /**********************
        Compute K5 coefficients
        */
        for (int i = 0; i < Ny_loc; i++) {
            shift_i = i*Nx_loc;
            for (int j = 0; j < Nx_loc; j++) {  
                lap_k4[shift_i+j] = F_parallel(K4, top, bottom, left, right, i, j, dx, dy, Ny_loc, Nx_loc, i + pos_i, j + pos_j, Ny, Nx);
                K5[shift_i+j] = tau*(lap_u[shift_i+j] + 0.5*lap_k1[shift_i+j]-1.5*lap_k3[shift_i+j] +
                                     2*lap_k4[shift_i+j]);
            }
        }
    
        // Declare array of K's
        double *K[5] = {K1, K2, K3, K4, K5};

        // Compute error in each subdomain, find max error w/ Reduce, bcast error
        E = max_in_arr(K, loc_arr_len);
        MPI_Reduce(&E, &E_global, 1, MPI_DOUBLE, MPI_MAX, root, comm);
        MPI_Bcast(&E_global, 1, MPI_DOUBLE, root, comm);
            
        // Update solution in each subdomain
        if(E_global < delta) {
            update_u_contiguous(local_u, K, loc_arr_len); 
            t += tau;
            
            // Gather u in each subdomain & write to file
            if (iproc == root) {
                collect_and_write_u(comm, local_u, m1, m2, loc_arr_len, t, x, y, Ny, Nx, root, filename, &f);
            } else {
                // Other processess do nothing
                MPI_Gather(local_u, loc_arr_len, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, root, comm);
            }
            
        }

        // Compute new value of the timestep
        tau = pow(delta/E_global, 0.2)*omega*tau;

        //MPI_Bcast(&last, 1, MPI_C_BOOL, root, comm);


        if (last) 
            break;
        
            
        #ifdef DEBUG


        if (iproc == 0) {
            std::cout << "E = " << E_global << std::endl;
            std::cout << "t = " << t << " tau = " << tau << " last = " << last << std::endl;
        }   
        //break;
    
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
    delete[] top;
    delete[] bottom;
    delete[] left;
    delete[] right;

}

void convolution_in_t(int Ny, int Nx, double *x, double *y, double **u, double t) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            u[i][j] = convolve(Ny, Nx, x, y, t, x[j], y[i]);
        }
    }
}

double convolve(int Nx, int Ny, double *x, double *y, double t, double m, double n) {
    double conv = 0;
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];

    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            conv += Gauss(m-x[i], n-y[j], t)*(signum(-sqrt(x[j]*x[j] + y[i]*y[i]) + 0.1)+1)*dx*dy;
        }
    }
    return conv;
}
