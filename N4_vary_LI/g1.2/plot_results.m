%%
clear;clc;
addpath('../../func_aux');
%% Original Image
% figure();hold all;
% thm.special = {};
% coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% 
% coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_anti', 'lambdaI', 'MAX(x)');
% coco_plot_bd(thm, 'run_zero', 'lambdaI', 'x');
% 
% box on;
% 
% title('g = 0.1 / N = 4');
% ylabel('R_1^*');
% exportgraphics(gcf, 'g0.1.eps', 'Resolution', 300);


%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');


[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~, x] = find_special_points('run_zero', 'HB');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


plot_curve_po('run_anti', colors(4,:), '--', 0.5);


xlim([-6, 1]);
ylim([-0.1, 1.11]);
box on;
xlabel('\lambda_i');
ylabel('R_1')
title('g = 1.2');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g1.2.eps', 'Resolution', 300);

