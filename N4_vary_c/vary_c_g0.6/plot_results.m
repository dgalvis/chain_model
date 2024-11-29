%%
clear;clc;
addpath('../../func_aux');
%% Original Image
% figure();hold all;
% thm = struct();
% thm.special = {'BP'};
% %thm.lspec = {'g-', 'LineWidth', 1};
% %coco_plot_bd(thm, 'run_sync_c', 'lambdaI', 'x');
% 
% 
% coco_plot_bd(thm, 'run_symm_c', 'lambdaI', 'x');
% coco_plot_bd(thm, 'run_anti_c', 'lambdaI', 'x');
% 
% box on;
% 
% title('g = 0.6  / N = 4');
% ylabel('R_1^*');

%%
figure();hold all;
colors = get_colors();


%plot_curve('run_sync_c', colors(1,:), '--', 'cc');
plot_curve('run_symm_c', colors(3,:), '--', 'cc');
plot_curve('run_anti_c', colors(4,:), '--', 'cc'); 
plot_curve('run_loop_A', colors(5,:), '--', 'cc'); 
plot_curve('run_loop_B', colors(5,:), '--', 'cc'); 
[~, ~, cc, x] = find_special_points('run_symm_c', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_anti_c', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end


xlim([-5, 5]);
%ylim([0.9, 1.01]);
box on;
xlabel('c');
ylabel('R_1')
title('g = 0.6, \lambda_i = -1.0');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, 'g0.6_vary_c.eps', 'Resolution', 300);
