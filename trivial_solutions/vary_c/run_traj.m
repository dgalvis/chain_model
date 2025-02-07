%% Initialise
clear; clc; close all;
addpath('../../functions'); % Add the 'functions' directory to the search path

%% Initial parameter values and Fsolve information


N = 4;
% Fixed parameters
LE = 1;                % Length of the element
c  = 5;                % Constant parameter
omega = 2*pi;          % Angular frequency

% Varied parameters
LI = 0.2;              % Initial value for the length parameter (will be varied)
g  = 2;              % Initial value for the gain parameter (will be varied)


p = [LE; LI; g; c; omega];

%% Trivial
rng(1);
R0 = rand(1, N);
T0 = [0, pi*rand(1, N-1)];
y0 = [R0.*cos(T0), R0.*sin(T0)];
%y0 = rand(1, 2*N);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, p, N));
[t,y] = ode15s(@(t,y)SL_vf(t, y, p, N), tspan, y0, options);

R = sqrt(y(:, 1:N).^2 + y(:, (N+1):end).^2);
figure();plot(R);pause(0.1);


rng(10);
R0 = rand(1, N);
T0 = [0, pi*rand(1, N-1)];
y0 = [R0.*cos(T0), R0.*sin(T0)];
%y0 = rand(1, 2*N);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, p, N));
[t,y] = ode15s(@(t,y)SL_vf(t, y, p, N), tspan, y0, options);

R = sqrt(y(:, 1:N).^2 + y(:, (N+1):end).^2);
figure();plot(R);pause(0.1);

