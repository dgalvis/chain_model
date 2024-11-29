%%
clear;clc;
addpath('../../func_aux');
%% Original Image
% figure();hold all;
% thm.special = {'HB'};
% coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% 
% coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_anti', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');
% thm.special = {};
% coco_plot_bd(thm, 'run_anti_HB_00', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_anti_HB_00', 'lambdaI', 'MIN(x)');
% 
% box on;
% 
% title('g = 0.1 / N = 4');
% ylabel('R_1^*');

%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');
plot_curve_po('run_anti_HB_00', colors(9,:), '--');

[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~, x] = find_special_points('run_anti', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

xlim([-1, 1]);
ylim([-0.1, 1.1]);
box on;
xlabel('lambdaI');
ylabel('R_1')
title('g = 0.6');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g0.6.eps', 'Resolution', 300);




