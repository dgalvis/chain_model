%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values

N = 4;
LE = 1;
LI = 1;
g  = 1.1;
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -6;
LI_lim = [LI_min, 6];
%% Synchronous solution continuation

y0 = [ones(1,N), zeros(1, N-1)]';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [3000,1000]);

coco(prob, 'run_sync', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%% Symmetric Branch

bd = coco_bd_read('run_sync');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,0]);
prob = ode_BP2ep(prob, '', 'run_sync', labs(2));
coco(prob, 'run_symm', [], 1, {'lambdaI', 'g'}, LI_lim);

%% Loop

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(2));
coco(prob, 'run_loop', [], 1, {'lambdaI', 'g'}, LI_lim);


%% Anti-symmetric wave

rng(1);
paux = p;
paux(2) = -2;
y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

y0 = y(end, :)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_anti', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, paux, 1, {'lambdaI', 'g'}, LI_lim);


%% Trivial Solution

y0 = zeros(1, 2*N)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_zero', @ode_isol2ep, @(y, p)SL_vf(0, y, p, N),@(y, p)SL_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);
