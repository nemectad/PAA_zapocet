#include "solvers.h"
#include "utilities.h"
#include <math.h>


double Laplace_u(double **u, int i, int j, double *x, double *y) {
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];
    return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx) + 
           (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy);
}

double F(double **u, int i, int j, double *x, double *y, int Nx, int Ny) {
    // Boundary condition
    if (i == 0 || j == 0 || i == Nx || j == Ny) {
        return 0;
    }
    // Laplace
    return Laplace_u(u, i, j, x, y);
}

void update_u(double **u, double ***K, int Nx, int Ny) {
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            u[i][j] = u[i][j] + (K[0][i][j] + 4*K[3][i][j] + K[4][i][j])/6;
        }
    }
}

void Merson(double **u, double *x, double *y, int Nx, int Ny, double dt, 
            double T, double delta) {
    double t, tau, eps, E, omega;
    bool last = false;
    tau = dt;
    std::ofstream f;
    std::string filename = "data.txt";
    omega = 0.8;
    t = 0.0;
    
    double **lap_u = new double*[Nx+1];
    double **lap_k1 = new double*[Nx+1];
    double **lap_k2 = new double*[Nx+1];
    double **lap_k3 = new double*[Nx+1];
    double **lap_k4 = new double*[Nx+1];
    double **lap_k5 = new double*[Nx+1];
    double ***K = new double**[5];
    // K matrix initialization - (5)*(Nx+1)*(Ny+1) array
    for (int i = 0; i < 5; i++) {
        K[i] = new double*[Nx+1];
        alloc(K[i], Nx, Ny);
    }

    alloc(lap_u, Nx, Ny);
    alloc(lap_k1, Nx, Ny);
    alloc(lap_k2, Nx, Ny);
    alloc(lap_k3, Nx, Ny);
    alloc(lap_k4, Nx, Ny);
    alloc(lap_k5, Nx, Ny);

    // Write initial data:
    f.open(filename, std::ios_base::app);
    f << Nx << "\t" << Ny << "\n";
    f.close();
    write_data(u, t, x, y, Nx, Ny, filename, &f);

    // Main loop
    while(1) {
        if (abs(T-t) < abs(tau)) {
            tau = T - t;
            last = true;
        }

        /*
        In this part we exploit the linearity of the laplacian. 
        We can therefore apply the function F (i.e. laplacian) and then add up all
        the terms correspondingly. 
        */ 
        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {  
                lap_u[i][j] = F(u, i, j, x, y, Nx, Ny);         
                K[0][i][j] = tau*lap_u[i][j];
            }
        }

        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {  
                lap_k1[i][j] = F(K[0], i, j, x, y, Nx, Ny);
                K[1][i][j] = tau*(lap_u[i][j] + lap_k1[i][j]/3);
            }
        }

        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {  
                lap_k2[i][j] = F(K[1], i, j, x, y, Nx, Ny);
                K[2][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + lap_k2[i][j])/6);
            }
        }

        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {  
                lap_k3[i][j] = F(K[2], i, j, x, y, Nx, Ny);
                K[3][i][j] = tau*(lap_u[i][j] + (lap_k1[i][j] + 3*lap_k3[i][j])/8);
            }
        }

        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {  
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

        //break;

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
    K = nullptr;
    
}