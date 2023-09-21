function x = icn(f, g, x, t, simParam)
% Iterative Crank-Nicolson method for time advancement of Nonlinear ODE
% dx/dt = f(x, t). For reference see: 
% Input
% Return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao and Thomas Bewley
% Institute    :    Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Sep. 21, 2023
if nargin < 5 || ~isfield(simParam, 'h')
    simParam.h = t + sign(t)*1e-10;
end

if nargin < 5 || ~isfield(simParam, 'T')
    simParam.T = simParam.h;
end

if nargin < 5 ||  ~isfield(simParam, 'method')
    simParam.method = 'Picard';
end

n = size(x, 1); h = simParam.h; T = simParam.T; method = simParam.method;
eps = 1e-3;
omega = .5;

while abs(t) < abs(T)
    acs = inf; % absolute change in solution
    ar = inf; % absolute residual
    xhat = x; % algorithm starts from the value of the previous time step
    while acs > eps && ar > eps 
        x_prev = xhat;
        if strcmp(method, 'Newton')
            % Newton method
            F = xhat - x - h/2*(f(xhat)+f(x));
            Fprime = eye(n) - h/2*g(xhat);
            xhat = xhat - omega*Fprime\F;
        elseif strcmp(method, 'Picard')
            % Picard 
            F = xhat - x - h/2*(f(xhat)+f(x));
            xstar = x + h/2*(f(xhat)+f(x));
            xhat = omega*xstar + (1-omega)*xhat;
        end
        ar = norm(F);
        acs = norm(xhat-x_prev);
    end
    t = t + h;
    x = xhat;
end
end
