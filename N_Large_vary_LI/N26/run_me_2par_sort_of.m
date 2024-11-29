%% Initialise
clear;clc;close all;
addpath('../../functions');
%% Initial parameter values
N = 26;

LE = 1;LI = 0.5;g  = 1;c  = 1;omega = 2*pi;
p = [LE; LI; g; c; omega];

LI_min = -0.25;
LI_lim = [LI_min, 0.5];

glist = 1:0.1:2.5;
figure();hold all;
for i = 1:length(glist)

    p(3) = glist(i);
    name_sync = ['run_sync_', num2str(i)];
    name_symm = ['run_symm_', num2str(i)];
    name_loop = ['run_loop_', num2str(i)];
    name_anti = ['run_anti_', num2str(i)];
    
    % Synchronous solution continuation
    y0 = [ones(1,N), zeros(1, N-1)]';
    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'vectorized', true);
    prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
    prob = coco_set(prob, 'cont',   'h0' , 0.01, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [1000,0]);
    coco(prob, name_sync, @ode_isol2ep, @(y, p)SL_polar_vf(0, y, p, N),@(y, p)SL_polar_jac(0, y, p, N), [], y0, ...
        {'lambdaE', 'lambdaI', 'g', 'c', 'omega'}, p, 1, {'lambdaI', 'g'}, LI_lim);
    
    % Symmetric Branch
    
    
    bd = coco_bd_read(name_sync);
    labs = coco_bd_labs(bd, 'BP');
    
    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'vectorized', true);
    prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
    prob = coco_set(prob, 'cont',   'h0' , 0.01, 'h_min', 1e-12, 'h_max', 0.05, 'NAdapt', 10, 'PtMX', [1000,0]);
    prob = ode_BP2ep(prob, '', name_sync, labs(end));
    coco(prob, name_symm, [], 1, {'lambdaI', 'g'}, LI_lim);

    % Loop
    
    bd = coco_bd_read(name_symm);
    labs = coco_bd_labs(bd, 'BP');
    
    
    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'vectorized', true);
    prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
    prob = coco_set(prob, 'cont',   'h0' , 0.001, 'h_min', 1e-12, 'h_max', 0.1, 'NAdapt', 10, 'PtMX', [0,300]);
    prob = ode_BP2ep(prob, '', name_symm, labs(end-1));
    coco(prob, name_loop, [], 1, {'lambdaI', 'g'}, LI_lim);

    %Anti
    bd = coco_bd_read(name_loop);
    labs = coco_bd_labs(bd, 'BP');
    
    prob = coco_prob();
    prob = coco_set(prob, 'ode', 'vectorized', true);
    prob = coco_set(prob, 'coll', 'NTST', 20, 'NCOL', 4);
    prob = coco_set(prob, 'cont',   'h0' , 0.05, 'h_min', 1e-6, 'h_max', 0.05, 'NAdapt', 10, 'PtMX', [1000,10]);
    prob = ode_BP2ep(prob, '', name_loop, labs(end-1));
    coco(prob, name_anti, [], 1, {'lambdaI', 'g'}, LI_lim);

    %Plot
    
    clf;hold all;
    thm.special = {'BP', 'HB'};
    coco_plot_bd(thm, name_sync, 'lambdaI', 'x');
    coco_plot_bd(thm, name_symm, 'lambdaI', 'x');
    coco_plot_bd(thm, name_anti, 'lambdaI', 'x');
    coco_plot_bd(thm, name_loop, 'lambdaI', 'x');
    pause(0.0001);

end

