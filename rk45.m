function [x, t] = rk45(f, x, t, simParam)
% rk45 time marching code, slightly modified from NR P292.
% Input:
% feval: the function handle, x'= f(x)
% x: the initial state
% t: the initial time
% h: the initla time step
% p: the time-dependent (or, state-depenedent) parameter for feval
% 
% Return:
% x: final state
% t: final time
if nargin < 4 || ~isfield(simParam, 'h')
    simParam.h = t+sign(t)*1e-10;
end

if nargin < 4 || ~isfield(simParam, 'T')
    simParam.T = simParam.h;
end

if nargin < 4 || ~isfield(simParam, 'epsOverT')
    simParam.epsOverT = 1e-2;
end

H = simParam.h; h = H/2; T = simParam.T;
while abs(t) < abs(T)
    f1 = feval(f, x); 
    f2 = feval(f, x + H * f1/2); 
    f3 = feval(f, x + H * f2/2);
    f4 = feval(f, x + H * f3);
    X = x + H * (f1/6 + (f2 + f3)/3 + f4/6);

    f1 = feval(f, x); 
    f2 = feval(f, x + h * f1/2);
    f3 = feval(f, x + h * f2/2);
    f4 = feval(f, x + h * f3);
    x = x + h * (f1/6 + (f2 + f3)/3 + f4/6);
    delta = norm(x - X, 1) / 15;

    x = (x * 16 - X) / 15; t = t + H;
    H = min(abs(H)*(abs(H)*simParam.epsOverT / delta)^(1/4), abs(simParam.T-t));
end
end
