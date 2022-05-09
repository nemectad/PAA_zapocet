#include "solvers.h"
#include "utilities.h"
#include <cmath>



double Laplace_u(double **u, int i, int j, double *x, double *y) {
    double dx = x[1]-x[0];
    double dy = y[1]-y[0];
    return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx) + 
           (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy);
}

double F(double **u, int i, int j, double *x, double *y, int Nx, int Ny) {
    // Boundary condition
    if (i == 0 || j == 0 || i == Nx-1 || j == Ny-1) {
        return 0;
    }
    // Laplace
    return Laplace_u(u, i, j, x, y);
}

void update_u(double **u, double ***K, int Nx, int Ny) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = u[i][j] + (K[0][i][j] + 4*K[3][i][j] + K[4][i][j])/6;
        }
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

void Merson_parallel(int argc, char* argv[], double **u, double *x, double *y, int Nx, int Ny, double dt, 
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
    /*
    f.open(filename, std::ios::out | std::ios_base::app);
    f << Nx << "\t" << Ny << "\n";
    f.close();
    write_data(u, t, x, y, Nx, Ny, filename, &f);

    */
    // Main loop
    while(1) {
        if (std::abs(T-t) < std::abs(tau)) {
            tau = T - t;
            last = true;
        }

        int nproc, iproc;
        
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

        
        // Root process
        int root = 0;
        int m1, m2;
        // Separate the main domain into m1*m2 subdomains
        m1 = 2;
        m2 = 2;

        int stride = (Nx/m1)*(Ny/m2);

        // Define custom MPI type => (Nx/m1)x(Ny/m2) 1D vector
        MPI_Datatype stype;

        MPI_Type_contiguous(stride, MPI_DOUBLE, &stype);
        MPI_Type_commit(&stype);

        // Broadcast the variables across all processes
        //int succ = MPI_Bcast(&Nx, 1, MPI_INT, root, MPI_COMM_WORLD);
        //std::cout << u[25][25] << std::endl;
        
        /*
        if (iproc != root) {
            std::cout << "m1 = " << m1 << std::endl;
            std::cout << "Nx = " << Nx << " in process " << iproc << std::endl;
            
        }
        */

        
        /*
        In this part we exploit the linearity of the laplacian. 
        We can therefore apply the function F (i.e. laplacian) and then add up all
        the terms correspondingly. 
        */ 

        // Sem by mel prijit broadcast/scatter hodnoty u, gather K, pak broadcast/scatter K[0] apod.

        // Declare local array for computation on the subdomain (has the same 
        // number of elements as the corresponding matrix)
        double local_K[stride];
        double local_lap[stride];
        //alloc(local_K, Nx/m1, Ny/m2);
        //alloc(local_lap, Nx/m1, Ny/m2);

        /*
        for (int i = (iproc%m1)*(Nx/m1); i < Nx/(m1-iproc%m1); i++) {
            for (int j = (iproc%m2)*(Nx/m2); j < Ny/(m2-iproc%m2); j++) {  
                local_lap[i*(Nx/2)+j] = F(u, i, j, x, y, Nx, Ny);         
                //lap_u[i][j] = F(u, i, j, x, y, Nx, Ny);         
                //K[0][i][j] = tau*lap_u[i][j];
                local_K[i*(Nx/2)+j] = tau*local_lap[i*(Nx/2)+j];
            }
        }
        */
        for (int i = 0; i < Nx/m1; i++) {
            int shift_i = i*(Nx/m1);
            int pos_i = (iproc/m1)*(Nx/m1);
            int pos_j = (iproc%m2)*(Ny/m2);
            for (int j = 0; j < Ny/m2; j++) {  
                local_lap[shift_i+j] = F(u, i + pos_i, j + pos_j, x, y, Nx, Ny);         
                //lap_u[i][j] = F(u, i, j, x, y, Nx, Ny);         
                //K[0][i][j] = tau*lap_u[i][j];
                local_K[shift_i+j] = tau*local_lap[shift_i+j];
            }
        }
        
        //std::cout << local_K[stride-1] << std::endl;

        if (iproc == root) {
            int size = m1*m2;

            //std::cout << stride << std::endl;

            double *buffer_K = new double[size*stride];
            double *buffer_lap = new double[size*stride];

            /*for (int i = 0; i<size; i++) {
                displacements[i] = i*stride;
                count_recv[i] = (Nx/m1)*(Ny/m2);
            }
            */
            
            int success = MPI_Gather(local_K, stride, MPI_DOUBLE, buffer_K, 1, stype, root, MPI_COMM_WORLD);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, buffer_lap, 1, stype, root, MPI_COMM_WORLD);
            //int success = MPI_Gatherv(local_K, 1, stype, buffer_K, count_recv, displacements, MPI_DOUBLE, root, MPI_COMM_WORLD);
            //std::cout << success << std::endl;

            
            // Collect the flattened buffers
            for (int k1 = 0; k1 < m1; k1++) {
                int shift_k1 = k1*(Nx/m1);
                for (int k2 = 0; k2 < m2; k2++) {
                    int shift_k2 = k2*(Ny/m2);
                    for (int i = 0; i < Nx/m1; i++) {
                        int shift_i = i*(Nx/m1);
                        for (int j = 0; j < Ny/m2; j++) {
                            lap_u[i+shift_k1][j+shift_k2] = buffer_lap[shift_i+j];
                            K[0][i+shift_k1][j+shift_k2] = buffer_K[shift_i+j];
                        }
                    }
                }
            }

            delete[] buffer_K;
            delete[] buffer_lap;

            //std::cout << lap_u[20][20] << std::endl;
        } else {
            int success = MPI_Gather(local_K, stride, MPI_DOUBLE, NULL, 0, stype, root, MPI_COMM_WORLD);
            MPI_Gather(local_lap, stride, MPI_DOUBLE, NULL, 0, stype, root, MPI_COMM_WORLD);
        }

        

        //dealloc(local_K, Nx/m1);
        //dealloc(local_lap, Nx/m1);
        MPI_Finalize();
        break;


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


        //dealloc(local_K, Nx/m1);
        //dealloc(local_lap, Nx/m1);

        MPI_Finalize();

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

    
}

void convolution_in_t(int Nx, int Ny, double *x, double *y, double **u, double t) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            /*if (i == 0 || j == 0 || i == Nx || j == Ny) {
                u[i][j] = 0;
            }
            else*/
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
            /*if (i == 0 || j == 0 || i == Nx || j == Ny) {
                conv = conv;
            }
            else
            */
            conv += Gauss(m-x[i], n-y[j], t)*(signum(-sqrt(x[i]*x[i] + y[j]*y[j]) + 0.1)+1)*dx*dy;
        }
    }
    return conv;
}
