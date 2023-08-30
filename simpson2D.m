function int_f = simpson2D(f, Nx, dx)
% Computes the approximation of double integral of f using simpson's method with
% unit spacing gridpoints Nx-by-Ny; notice that Nx and Ny are the number of 
% gridpoints for x and y axis respectively, thus must be odd.
%
% Input:
% f: vector
% Nx, Ny, dx and dy: scalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao
% Institute    :    Mechanical and Aerospace Engineering, UC San Diego
% Date         :    July. 04, 2023

if isscalar(Nx)
    Ny = Nx;
elseif isvector(Nx)
    Ny = Nx(2);
    Nx = Nx(1);
else
    error("Input 'Nx' should either be scalar or vector. ");
end

if isscalar(dx)
    dy = dx;
elseif isvector(dx)
    dy = dx(2);
    dx = dx(1);
else
    error("Input 'dx' should either be scalar or vector. ");
end

if mod(Nx, 2) ~= 1
    error('X-axis should have odd number of gridpoints. ');
elseif mod(Ny, 2) ~= 1
    error('Y-axis should have odd number of gridpoints. ');
end

if ~isequal(size(f), [Nx*Ny, 1])
    error('size of f is not equal to [Nx*Ny, 1]. ')
end
xx = [1, repmat([4, 2], 1, (Nx-3)/2), 4, 1];
yy = [1, repmat([4, 2], 1, (Ny-3)/2), 4, 1];
S = kron(xx, yy);
int_f = S * f * dx * dy / 9;
end
