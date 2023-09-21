function [x, t] = rk45(f, x, t, simParam)
% This script implements the classic four-stage, fifth-order accurate Rungze 
% Kutta scheme for time advancment of ODE system, dx/dt = f(x,t), with a slight
% modification of Prof. Thomas Bewley's Numerical Renaissance P292. This
% script is suitable to march the state both forward and backward in time.
% For backward marching, simply set t to be negative real number.
% Input:
%       f         :  function handle, x'= f(x)
%       x         :  initial state
%       t         :  initial time
%       simParam.h:  time step for marching
%       simParam.p:  time(or, state)-dependent parameter of the RHS of ODE
% 
% Return:
%       x         :  final state
%       t         :  final time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       :    Muhan Zhao and Thomas Bewley
% Institute    :    Mechanical and Aerospace Engineering, UC San Diego
% Date         :    Created Feb., 2023; Modified Sep. 21, 2023


%%%%%%%%%%
% TODO record x for all time instants b/w t and T
%%%%%%%%%%

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
    % calculate X with one RK4 full timestep H
    f1 = feval(f, x); 
    f2 = feval(f, x + H * f1/2); 
    f3 = feval(f, x + H * f2/2);
    f4 = feval(f, x + H * f3);
    X = x + H * (f1/6 + (f2 + f3)/3 + f4/6);

    % calculate X with two RK4 sub timesteps h=H/2
    f1 = feval(f, x); 
    f2 = feval(f, x + h * f1/2);
    f3 = feval(f, x + h * f2/2);
    f4 = feval(f, x + h * f3);
    x = x + h * (f1/6 + (f2 + f3)/3 + f4/6);
    f1 = feval(f, x); 
    f2 = feval(f, x + h * f1/2);
    f3 = feval(f, x + h * f2/2);
    f4 = feval(f, x + h * f3);
    x = x + h * (f1/6 + (f2 + f3)/3 + f4/6);

    % estimate of error
    delta = norm(x - X, 1) / 15;
    % determine the new x based on adaptive rk45
    x = (x * 16 - X) / 15; t = t + H;
    H = min(abs(H)*(abs(H)*simParam.epsOverT / delta)^(1/4), abs(simParam.T-t));
end
end
