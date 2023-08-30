function [NABx, NABy] = twoD_nabla_uniform(Nx, Ny, delta)
% Compute the nabla operator for x and y dimension, using second order
% central difference in the interior of the domain, and second order
% forward/backward difference at the boundary of the domain.
% 
% This script aims to accelerate the computation via the avoidance of the 
% double for loop. Rather than that, the script uses kronecker matrix product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao
% Institute    :    Flow & Control Lab, Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Jul. 15, 2023

% 2nd order central differencing scheme
% Nabla operator for x-axis
ex = ones(Nx, 1);
Lx = spdiags([-ex/2, zeros(Nx, 1), ex/2], [-1, 0, 1], Nx, Nx);

% define boundary condition for x-axis

% 2nd forward difference
temp = zeros(1, Nx); temp(1) = -3/2; temp(2) = 2; temp(3) = -1/2;
Lx(1, :) = temp;
% 2nd backward difference
temp = zeros(1, Nx); temp(end) = 3/2; temp(end-1) = -2; temp(end-2) = 1/2;
Lx(end, :) = temp;

% assemble NABx operator 
NABx = kron(eye(Ny), Lx); NABx = NABx/delta;

% Nabla operator for y-axis
ey = ones(Nx, 1);
Ly = spdiags([-ey/2, ey/2], [-1, 1], Nx, Ny);

NABy = kron(Ly, eye(Nx));

% define boundary condition for y-axis

    % 2nd order forward difference
temp = cat(2, -3/2*eye(Ny), 2*eye(Ny), -1/2*eye(Ny));
NABy(1:Ny, 1:3*Ny) = temp;

% 2nd order backward difference
temp = cat(2, 1/2*eye(Ny), -2*eye(Ny), 3/2*eye(Ny));
NABy(end-Ny+1:end, end-3*Ny+1:end) = temp;

NABy = NABy/delta;

end
