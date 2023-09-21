function [NABx, NABy] = twoD_nabla_uniform(Nx, Ny, delta)
% Compute the nabla operator for x and y axis assuming a uniform spacing.
% The first order derivative of the interior of the domain is approximated
% by the second order central differencing scheme, while the second order
% forward/backward difference is applied at the boundary of the domain.
% 
% Both nabla operators are "row-wise", namely NAB_(:,j) denotes the value
% of discretized x-axis associated with y_j. This convention is compatible
% with the common understanding of the theoretical development in numerical 
% analysis.
% 
% User might define the boundary condition at their own interest with a
% slight modification in this script.
% 
% Input:
%       Nx      : The number of discretized points in x
%       Ny      : The number of discretized points in y
%       delta   : The step size of each axis
% Example:
% [NABx, NABy] = twoD_nabla_uniform(64, 64, .01);
% 
% This script aims to accelerate the computation by the substituion of the 
% double for loop with the kronecker matrix product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao
% Institute    :    Flow & Control Lab, Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Jul. 15, 2023

% 2nd order central differencing scheme
% Nabla operator for x-axis
ex = ones(Nx, 1);
Lx = spdiags([-ex/2, zeros(Nx, 1), ex/2], [-1, 0, 1], Nx, Nx);

% define boundary condition for x-axis

% 2nd forward difference for west boundary
temp = zeros(1, Nx); temp(1) = -3/2; temp(2) = 2; temp(3) = -1/2;
Lx(1, :) = temp;
% 2nd backward difference for east boundary
temp = zeros(1, Nx); temp(end) = 3/2; temp(end-1) = -2; temp(end-2) = 1/2;
Lx(end, :) = temp;

% assemble NABx operator 
NABx = kron(eye(Ny), Lx); NABx = NABx/delta;

% Nabla operator for y-axis
ey = ones(Nx, 1);
Ly = spdiags([-ey/2, ey/2], [-1, 1], Nx, Ny);

NABy = kron(Ly, eye(Nx));

% define boundary condition for y-axis

% 2nd order forward difference for south boundary
temp = cat(2, -3/2*eye(Ny), 2*eye(Ny), -1/2*eye(Ny));
NABy(1:Ny, 1:3*Ny) = temp;

% 2nd order backward difference for north boundary
temp = cat(2, 1/2*eye(Ny), -2*eye(Ny), 3/2*eye(Ny));
NABy(end-Ny+1:end, end-3*Ny+1:end) = temp;

NABy = NABy/delta;

end
