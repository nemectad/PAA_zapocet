function [t, X, Y, u] = sol(f, Ny, Nx)
    % Time
    t = textscan(fgetl(f), "%f"){1};
    % Meshgrid allocation
    X = zeros(Ny, Nx);
    Y = zeros(Ny, Nx);
    % Solution allocation
    u = zeros(Ny, Nx);
    for i = 1:(Ny)
        for j = 1:(Nx)
            line = textscan(fgetl(f), "%f %f %f");
            % Generate meshgrid for surface plot
            Y(i, j) = line{1};
            X(i, j) = line{2};
            % Obtain solution at each individual point [X, Y]
            u(i, j) = line{3};
        end
    end
    
end