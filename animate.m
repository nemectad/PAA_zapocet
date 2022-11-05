filename = input("Enter the name of the data file: ", "s");

f = fopen(filename);
line = textscan(fgetl(f), "%f %f");
Ny = line{1};
Nx = line{2};

[t, X, Y, u] = sol(f, Ny, Nx);

figure
mesh(X, Y, u);

pause

% Plot animation of the solution
while (! feof (f) )
    [t, X, Y, u] = sol(f, Ny, Nx);
    mesh(X, Y, u);
    pause(0.1)
end