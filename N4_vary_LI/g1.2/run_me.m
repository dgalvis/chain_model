%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values

N = 4;
LE = 1;
LI = 1;
g  = 1.2;
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -6;
LI_lim = [LI_min, 2];
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


%% Trivial Solution

y0 = zeros(1, 2*N)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_zero', @ode_isol2ep, @(y, p)SL_vf(0, y, p, N),@(y, p)SL_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%%
bd = coco_bd_read('run_zero');
labs = coco_bd_labs(bd, 'HB');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [0,1000]);
prob = ode_HB2po(prob, '', 'run_zero', labs(2));
coco(prob, 'run_anti', [], 1, {'lambdaI', 'g'}, LI_lim);

