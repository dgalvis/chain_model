%%
clear;close all;clc;
addpath('../func_aux');
% %%
% figure();hold all;
% thm = struct('special', {{'BP'}});
% coco_plot_bd(thm, 'run_anti_po', 'lambdaI', 'MAX(x)');
% thm = struct('special', {{'BP', 'HB'}});
% coco_plot_bd(thm, 'run_sync', 'lambdaI', 'x');
% thm = struct('special', {{'BP'}});
% coco_plot_bd(thm, 'run_symm', 'lambdaI', 'x');
% thm = struct('special', {{'BP'}});
% coco_plot_bd(thm, 'run_loop', 'lambdaI', 'x');
% 
% 
% title('N = 5, g = 0.1')
% ylabel('R_1^*')
% 
% box on;
% exportgraphics(gcf, 'N5_gsmall.eps', 'Resolution', 300);






%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
%plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_loop', colors(5,:), '--');


[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_loop', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


plot_curve_po('run_anti_po', colors(4,:), '--', 0.5);

axis tight
xlim([-1, 1]);
ylim([0.95, 1.01]);
box on;
xlabel('\lambda_{int}');
ylabel('R_1')
title('g = 0.1/N=5');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, 'g0.1.eps', 'Resolution', 300);

