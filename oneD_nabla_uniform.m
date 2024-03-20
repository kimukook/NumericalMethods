function Lx = oneD_nabla_uniform(Nx, delta, orderAccuracy)
% Compute the nabla operator for x axis assuming a uniform spacing discretization.
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
%       delta   : The step size of x-axis
% Example:
% NABx = twoD_nabla_uniform(64, 64, .01);
% 
% This script aims to accelerate the computation by the substituion of the 
% double for loop with the kronecker matrix product.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao
% Institute    :    Flow & Control Lab, Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Jul. 15, 2023

% 2nd-order central differencing scheme
% Nabla operator for x-axis
if nargin < 3
    orderAccuracy = 2;
end
switch orderAccuracy
    case 2
        ex = ones(Nx, 1); % unit vector with all entries being 1
        Lx = spdiags([-ex/2, zeros(Nx, 1), ex/2], [-1, 0, 1], Nx, Nx);
        
        % define gradient @ x-axis boundary 
        
    case 4
        ex = ones(Nx, 1);
        Lx = spdiags([ex/12, -2/3*ex, zeros(Nx, 1), 2/3*ex, -ex/12], [-2, -1, 0, 1, 2], Nx, Nx);
end
% 2nd forward difference for west boundary
temp = zeros(1, Nx); temp(1) = -3/2; temp(2) = 2; temp(3) = -1/2;
Lx(1, :) = temp;
% 2nd backward difference for east boundary
temp = zeros(1, Nx); temp(end) = 3/2; temp(end-1) = -2; temp(end-2) = 1/2;
Lx(end, :) = temp;
Lx = Lx/delta;
end
