%% Initialise
clear;clc;close all;
addpath('../functions');
%% Initial parameter values

N = 5;
LE = 1;
LI = 1;
g  = 0.1;
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -5;
LI_lim = [LI_min, 1];

%% Synchronous solution continuation

y0 = [ones(1,N), zeros(1, N-1)]';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_sync', @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);


%% Symmetric Branch

bd = coco_bd_read('run_sync');
labs = coco_bd_labs(bd, 'BP');

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,0]);
prob = ode_BP2ep(prob, '', 'run_sync', labs(end));
coco(prob, 'run_symm', [], 1, {'lambdaI', 'g'}, LI_lim);


%% Loop

bd = coco_bd_read('run_symm');
labs = coco_bd_labs(bd, 'BP');


prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [0,1000]);
prob = ode_BP2ep(prob, '', 'run_symm', labs(1));
coco(prob, 'run_loop', [], 1, {'lambdaI', 'g'}, LI_lim);


%% Antisymmetric branch

rng(1);
paux = p;
paux(2) = -0.1;

y0 = rand(1,2*N);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_vf(t,y,paux, N), tspan, y0, options);

pols = y(end, :)';
R = sqrt(pols(1:N).^2 + pols(N+1:end).^2);
T = atan2(pols(N+1:end), pols(1:N));

omega_hat = omega + c*paux(2)*(1-R(1)^2) + g * R(2)*sin(T(2)-T(1))/R(1);

tmax = 2*pi/omega_hat;
t = 0:0.01:tmax;
y = [R.*cos(omega_hat*t + T);R.*sin(omega_hat*t + T)]';

%%
prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,2000]);
coll_func = {@(y,p)SL_vf(0,y,p,N), @(y,p)SL_jac(0,y,p,N), []};
coll_args = [coll_func, {t, y, {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, paux}];
prob = ode_isol2po(prob, '', coll_args{:});
coco(prob, 'run_anti_po', [], 1, {'lambdaI', 'g', 'po.period   '}, LI_lim);


%% Trivial Solution

y0 = zeros(1, 2*N)';

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.01, 'NAdapt', 10, 'PtMX', [1000,1000]);
coco(prob, 'run_zero', @ode_isol2ep, @(y, p)SL_vf(0, y, p, N),@(y, p)SL_jac(0, y, p, N), [], y0, ...
    {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);






