%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values
N = 16;

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
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.006, 'NAdapt', 10, 'PtMX', [10000,0]);
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



%% Anti-symmetric wave

bd = coco_bd_read('run_loop');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [2000,2000]);
prob = ode_BP2ep(prob, '', 'run_loop', labs(1));
coco(prob, 'run_anti', [], 1, {'lambdaI', 'g'}, LI_lim);



%% Loop 2

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.01, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,1950]);
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

bd = coco_bd_read('run_anti');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.01, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [3000,3000]);
prob = ode_BP2ep(prob, '', 'run_anti', labs(5));
coco(prob, 'run_loop_03', [], 1, {'lambdaI', 'g'}, LI_lim);

%%
% %% Symmetric wave HB
% 
% bd = coco_bd_read('run_symm');
% labs = coco_bd_labs(bd, 'HB');
% 
% prob = coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true);
% prob = coco_set(prob, 'coll', 'NTST',4, 'NCOL', 4);
% prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 1, 'NAdapt', 10, 'PtMX', [0, 5000]);
% prob = ode_HB2po(prob, '', 'run_symm', labs(2));
% coco(prob, 'run_symm_HB_00', [], 1, {'lambdaI', 'po.period'}, LI_lim);
% 
% 
% %% Anti-symmetric wave HB
% 
% bd = coco_bd_read('run_anti');
% labs = coco_bd_labs(bd, 'HB');
% 
% prob = coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true);
% prob = coco_set(prob, 'coll', 'NTST', 10, 'NCOL', 4);
% prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max',1, 'NAdapt', 10, 'PtMX', [0, 5000]);
% prob = ode_HB2po(prob, '', 'run_anti', labs(2));
% coco(prob, 'run_anti_HB_00', [], 1, {'lambdaI', 'po.period'}, LI_lim);
% 
% %% Loop HB
% bd = coco_bd_read('run_loop_03');
% labs = coco_bd_labs(bd, 'HB');
% 
% prob = coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true);
% prob = coco_set(prob, 'coll', 'NTST', 10, 'NCOL', 4);
% prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 1, 'NAdapt', 10, 'PtMX', [0, 5000]);
% prob = ode_HB2po(prob, '', 'run_loop_03', labs(2));
% coco(prob, 'run_loop_HB_00', [], 1, {'lambdaI', 'po.period'}, LI_lim);

%%
%
% figure();hold all;
% thm.special = {'BP', 'HB'};
% coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_anti', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop_02', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop_03', 'lambdaI', 'x');
% 
% % thm.special = {'BP'};
% % coco_plot_bd(thm, 'run_anti_HB_00', 'lambdaI', 'MAX(x)');
% % coco_plot_bd(thm, 'run_anti_HB_00', 'lambdaI', 'MIN(x)');
% % coco_plot_bd(thm, 'run_symm_HB_00', 'lambdaI', 'MAX(x)');
% % coco_plot_bd(thm, 'run_symm_HB_00', 'lambdaI', 'MIN(x)');
% % coco_plot_bd(thm, 'run_loop_HB_00', 'lambdaI', 'MAX(x)');
% % coco_plot_bd(thm, 'run_loop_HB_00', 'lambdaI', 'MIN(x)');
% ylim([-0.1, 1.1]);
% xlim([-0.5, 0.3]);

