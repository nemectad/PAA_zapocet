function [t, X, Y, u] = sol(f, Nx, Ny)
    % Time
    t = textscan(fgetl(f), "%f"){1};
    % Meshgrid allocation
    X = zeros(Nx, Ny);
    Y = zeros(Nx, Ny);
    % Solution allocation
    u = zeros(Nx, Ny);
    for i = 1:(Nx)
        for j = 1:(Ny)
            line = textscan(fgetl(f), "%f %f %f");
            % Generate meshgrid for surface plot
            X(i, j) = line{1};
            Y(i, j) = line{2};
            % Obtain solution at each individual point [X, Y]
            u(i, j) = line{3};
        end
    end
    
end