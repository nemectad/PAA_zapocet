function [t, X, Y, u] = sol(f, Nx, Ny)
    % Time
    t = textscan(fgetl(f), "%f"){1};
    % Meshgrid allocation
    X = zeros(Nx+1, Ny+1);
    Y = zeros(Nx+1, Ny+1);
    % Solution allocation
    u = zeros(Nx+1, Ny+1);
    for i = 1:Nx+1
        for j = 1:Ny+1
            line = textscan(fgetl(f), "%f %f %f");
            % Generate meshgrid for surface plot
            X(i, j) = line{1};
            Y(i, j) = line{2};
            % Obtain solution at each individual point [X, Y]
            u(i, j) = line{3};
        end
    end
    
end