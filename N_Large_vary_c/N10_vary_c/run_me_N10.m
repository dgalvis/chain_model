%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values
N = 10;

LE = 1;LI = 1;g  = 0.9;c  = 1;omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -0.2;
LI_lim = [LI_min, 1];
c_lim = [-2, 2];
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
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [1000,0]);
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
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [2000,2000]);
prob = ode_BP2ep(prob, '', 'run_loop', labs(end-1));
coco(prob, 'run_anti', [], 1, {'lambdaI', 'g'}, LI_lim);


%%
bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'EP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.0001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [10000,10000]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(1));
coco(prob, 'run_symm_c', [], 1, {'c', 'lambdaI'}, c_lim);

%%
bd = coco_bd_read('run_anti');
labs = coco_bd_labs(bd, 'EP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.0001, 'h_min', 1e-12, 'h_max', 0.005, 'NAdapt', 10, 'PtMX', [10000,10000]);
prob = ode_BP2ep(prob, '', 'run_anti', labs(1));
coco(prob, 'run_anti_c', [], 1, {'c', 'lambdaI'}, c_lim);

%%
bd = coco_bd_read('run_sync');
labs = coco_bd_labs(bd, 'EP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [1000,1000]);
prob = ode_BP2ep(prob, '', 'run_sync', labs(1));
coco(prob, 'run_sync_c', [], 1, {'c', 'lambdaI'}, c_lim);


%%
bd = coco_bd_read('run_symm_c');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,1500]);
prob = ode_BP2ep(prob, '', 'run_symm_c', labs(end));
coco(prob, 'run_loop_c_2', [], 1, {'c', 'lambdaI'}, c_lim);


%%
bd = coco_bd_read('run_symm_c');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0, 1500]);
prob = ode_BP2ep(prob, '', 'run_symm_c', labs(end-1));
coco(prob, 'run_loop_c_1', [], 1, {'c', 'lambdaI'}, c_lim);

%%
bd = coco_bd_read('run_symm_c');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,1500]);
prob = ode_BP2ep(prob, '', 'run_symm_c', labs(end-2));
coco(prob, 'run_loop_c_3', [], 1, {'c', 'lambdaI'}, c_lim);


%%
bd = coco_bd_read('run_symm_c');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0, 1500]);
prob = ode_BP2ep(prob, '', 'run_symm_c', labs(end-3));
coco(prob, 'run_loop_c_4', [], 1, {'c', 'lambdaI'}, c_lim);