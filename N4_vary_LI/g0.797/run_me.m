%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values

N = 4;
LE = 1;
LI = -2;
g  = 0.797; %0.797, 0.798
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -5;
LI_lim = [LI_min, 100];

%%
y0 = [ones(1,N), zeros(1, N-1)]';
prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_sync', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%% Find the Anti-symmetric wave

% Hi Kyle, hopefully your seed is the same as mine
rng(1);
y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, p, N));
[t,y] = ode15s(@(t,y)SL_polar_vf(t, y, p, N), tspan, y0, options);

y0 = y(end, :)'; 

err1 = abs(y0(1) - y0(4));
err2 = abs(y0(2) - y0(3));
err3 = abs(y0(end) - pi);
err4 = abs(y0(end-1) - y0(end-2) - pi);
disp('if all errs are 0, then asym wave: ')
errs = [err1, err2, err3, err4]


%% 1-D continuation

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [10000,10000]);
coco(prob, 'run_anti', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);


%% Loop
bd = coco_bd_read('run_anti');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [0,150]);
prob = ode_BP2ep(prob, '', 'run_anti', labs(1));
coco(prob, 'run_loop', [], 1, {'lambdaI', 'g'}, LI_lim);




%% Find the symmetric wave

paux = p;
paux(2) = -0.5;
% Hi Kyle, hopefully your seed is the same as mine
rng(1);
y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

y0 = y(end, :)'; 

err1 = abs(y0(1) - y0(4));
err2 = abs(y0(2) - y0(3));
err3 = abs(y0(end));
err4 = abs(y0(end-1) - y0(end-2));
disp('if all errs are 0, then asym wave: ')
errs = [err1, err2, err3, err4]

%% 1-D continuation

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [10000,1000]);
coco(prob, 'run_symm', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);


%% Loop
% bd = coco_bd_read('run_symm');
% labs = coco_bd_labs(bd, 'BP');
% 
% prob = coco_prob();
% prob = coco_set(prob, 'ode', 'vectorized', true);
% prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [0,100]);
% prob = ode_BP2ep(prob, '', 'run_symm', labs(1));
% coco(prob, 'run_loop_B', [], 1, {'lambdaI', 'g'}, LI_lim);

%% Trivial Solution

y0 = zeros(1, 2*N)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_zero', @ode_isol2ep, @(y, p)SL_vf(0, y, p, N),@(y, p)SL_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);

%%
figure();hold all;
thm = struct('special', {{'SN', 'HB'}});
coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
%coco_plot_bd(thm, 'run_loop_B', 'lambdaI', 'x');
coco_plot_bd(thm, 'run_anti', 'lambdaI', 'x');
coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');