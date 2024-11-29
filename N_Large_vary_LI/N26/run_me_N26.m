%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values
N = 26;

LE = 1;LI = 1;g  = 2.5;c  = 1;omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -1;
LI_lim = [LI_min, 1];
%% Synchronous solution continuation

y0 = [ones(1,N), zeros(1, N-1)]';
prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [5000,5000]);
coco(prob, 'run_sync', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%% Symmetric Branch


bd = coco_bd_read('run_sync');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [2000,0]);
prob = ode_BP2ep(prob, '', 'run_sync', labs(end));
coco(prob, 'run_symm', [], 1, {'lambdaI', 'g'}, LI_lim);

%% Loop

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');


prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,3231]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(end-1));
coco(prob, 'run_loop', [], 1, {'lambdaI', 'g'}, LI_lim);


%% Loop 2

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.01, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,2500]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(end-2));
coco(prob, 'run_loop_02', [], 1, {'lambdaI', 'g'}, LI_lim);


%% Trivial Solution

y0 = zeros(1, 2*N)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.001, 'NAdapt', 10, 'PtMX', [10000,10000]);
coco(prob, 'run_zero', @ode_isol2ep, @(y, p)SL_vf(0, y, p, N),@(y, p)SL_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%% Loop 3

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,2500]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(end-3));
coco(prob, 'run_loop_03', [], 1, {'lambdaI', 'g'}, LI_lim);



%%
bd = coco_bd_read('run_loop');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.005, 'NAdapt', 10, 'PtMX', [20000,20000]);
prob = ode_BP2ep(prob, '', 'run_loop', labs(end-1));
coco(prob, 'run_anti', [], 1, {'lambdaI', 'g'}, LI_lim);

%%
%
% figure();hold all;
% thm.special = {'BP', 'HB'};
% coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_anti', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop_02', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop_03', 'lambdaI', 'x');


