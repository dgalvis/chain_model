%% Initialise
clear;clc;close all;
addpath('../functions');

%% Parameters
N = 5;
LE = 1;LI = 1;g  = 0.1;c  = 1;omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -5;
LI_lim = [LI_min, 5];


%% Sync solution 
omega_hat = omega;
tmax = 2*pi/omega_hat;
t = 0:0.01:tmax;


R = ones(N,1);
T = zeros(N,1);


y = [R.*cos(omega_hat*t + T);R.*sin(omega_hat*t + T)]';

%% Sync
prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [2000,2000]);
coll_func = {@(y,p)SL_vf(0,y,p,N), @(y,p)SL_jac(0,y,p,N), []};
coll_args = [coll_func, {t', y, {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p}];
prob = ode_isol2po(prob, '', coll_args{:});
coco(prob, 'run_po_sync', [], 1, {'lambdaI','po.period'}, LI_lim);

%% Symm

bd = coco_bd_read('run_po_sync');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [2000,0]);


prob = ode_BP2po(prob, '', 'run_po_sync', labs(end));
coco(prob, 'run_po_symm', [], 1, {'lambdaI','po.period'}, LI_lim);

%% Loop

bd = coco_bd_read('run_po_symm');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [0,525]);


prob = ode_BP2po(prob, '', 'run_po_symm', labs(1));
coco(prob, 'run_po_loop', [], 1, {'lambdaI','po.period'}, LI_lim);

%% Anti

bd = coco_bd_read('run_po_loop');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [2000,2000]);


prob = ode_BP2po(prob, '', 'run_po_loop', labs(1));
coco(prob, 'run_po_anti', [], 1, {'lambdaI','po.period'}, LI_lim);

%%
figure();hold all;
thm.special = {'BP'};
coco_plot_bd(thm, 'run_po_sync', 'lambdaI', 'MAX(x)');
coco_plot_bd(thm, 'run_po_symm', 'lambdaI', 'MAX(x)');
coco_plot_bd(thm, 'run_po_loop', 'lambdaI', 'MAX(x)');
coco_plot_bd(thm, 'run_po_anti', 'lambdaI', 'MAX(x)');
%ylim([-0.1, 1.1]);