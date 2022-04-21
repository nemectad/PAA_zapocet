f = fopen("data.txt");
line = textscan(fgetl(f), "%f %f");
Nx = line{1};
Ny = line{2};

[t, X, Y, u] = sol(f, Nx, Ny);

figure
mesh(X, Y, u);

pause

% Plot animation of the solution
while (! feof (f) )
    [t, X, Y, u] = sol(f, Nx, Ny);
    mesh(X, Y, u);
    pause(0.01)
end