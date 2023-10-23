function [ind, NABn] = twoD_uniform_Neumann_bc(Ns)
% Define the linear operator for the boundary conditions used in CAPOPS project 
% Input:
% Output:
% NABzf: zero flux boundary condition used in rho formulation
% NABd: (homogeneous) Dirichlet boundary condition
% 
% Find the index for the grid points at the south, north, west and east
% boundaries.
% 
% Seperate x-bounday and y-boundary, as if they have different
% number of gridpoints, dimensions of rows being concatenated is not
% consistent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao
% Institute    :    Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Aug., modified Oct., 2023.
if isvector(Ns)
    Nx = Ns(1); Ny = Ns(2);
elseif isscalar(Ns)
    Nx = Ns; Ny = Ns;
end
south = 1 : Nx;
west = 1 + Nx*(0:Ny-1); 
% west = 1 + Nx*(1:Ny-2); % remove two corner points
east = Nx*(1:Ny); 
% east = Nx*(2:Ny-1); % remove two corner points
north = 1+Nx*(Ny-1) : Nx*Ny;

% Notice that four corner points will be defined twice: 
% for bottom two points, flux passing through south boundary will be removed.
% for top two points, flux passing through west and east boundaries will be
% removed.
ind = cat(2, south, reshape(cat(1, west, east), 1, 2*Ny), north);

NABn = zeros(2*(Nx+Ny), Nx*Ny);
for k = 1 : size(ind, 2)
    c = ind(k);
    n = ind(k) + Nx;
    w = ind(k) - 1;
    e = ind(k) + 1;
    s = ind(k) - Nx;

    if c <= Nx
        % south boundary, forward diff.
        NABn(k, c) = -1;
        NABn(k, n) =  1;
    
    elseif mod(c, Nx) == 1 % && c > Nx && c < 1+Nx*(Ny-1)
        % west boundary, forward diff.
        NABn(k, c) = -1;
        NABn(k, e) =  1;

    elseif mod(ind(k), Nx) == 0 % && c > Nx && c < 1+Nx*(Ny-1)
        % east bounday, backward diff.
        NABn(k, c) =  1;
        NABn(k, w) = -1;
        
    elseif ind(k) >= 1+Nx*(Ny-1)
        % north boundary, backward diff.
        NABn(k, c) =  1;
        NABn(k, s) = -1;

    end
end

end % function twoD_uniform_zeroFlux
