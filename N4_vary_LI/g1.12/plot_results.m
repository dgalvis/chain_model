%%
clear;clc;
addpath('../../func_aux');
%% Original Image
% figure();hold all;
% thm.special = {};
% %coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% 
% %coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
% %coco_plot_bd(thm, 'run_anti', 'lambdaI', 'x');
% %coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');
% 
% coco_plot_bd(thm, 'run_loop_HB_01', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_loop_HB_01', 'lambdaI', 'MIN(x)');
% 
% coco_plot_bd(thm, 'run_loop_HB_02', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_loop_HB_02', 'lambdaI', 'MIN(x)');
% 
% coco_plot_bd(thm, 'run_loop_HB_03', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_loop_HB_03', 'lambdaI', 'MIN(x)');
% 
% coco_plot_bd(thm, 'run_loop_HB_04', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_loop_HB_04', 'lambdaI', 'MIN(x)');
% 
% 
% box on;
% 
% ylabel('R_1^*');
%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');



% Branch points
[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(3), x(3), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(4), x(4), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_zero', 'HB');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


% Hopf bifurcations
[LI, ~, ~,x] = find_special_points('run_anti', 'HB');
scatter(LI(3), x(3), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(4), x(4), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

[LI, ~, ~,x] = find_special_points('run_loop', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(3), x(3), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(4), x(4), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

% Saddle-nodes
[LI, ~, ~,x] = find_special_points('run_loop', 'SN');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));

xlim([-6, 1]);
ylim([-0.1, 1.11]);
box on;
xlabel('lambdaI');
ylabel('R_1')
title('g = 1.12');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g1.12.eps', 'Resolution', 300);

%%
figure();hold all;
colors = get_colors();
plot_curve('run_loop', colors(5,:), '--');
plot_curve_po('run_loop_HB_01', colors(9,:), '--');
plot_curve_po('run_loop_HB_02', colors(9,:), '--');
plot_curve_po('run_loop_HB_03', colors(9,:), '--');
plot_curve_po('run_loop_HB_04', colors(9,:), '--');

% Hopf bifurcations

[LI, ~, ~,x] = find_special_points('run_loop', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(3), x(3), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(4), x(4), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

% Saddle-nodes
[LI, ~, ~, x] = find_special_points('run_loop', 'SN');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));

%xlim([-6, 1]);
%ylim([-0.1, 1.11]);
box on;
xlabel('lambdaI');
ylabel('R_1')
title('g = 1.12, Loop and Hopfs');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g1.12_Hopfs.eps', 'Resolution', 300);


