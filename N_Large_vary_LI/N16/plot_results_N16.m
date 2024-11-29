%%
clear;clc;
addpath('../../func_aux');

%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');
plot_curve('run_loop_02', colors(5,:), '--');
%plot_curve('run_loop_03', colors(5,:), '--');
%plot_curve_po('run_anti_HB_00', colors(9,:), '--'); 
%plot_curve_po('run_symm_HB_00', colors(9,:), '--'); 
%plot_curve_po('run_loop_HB_00', colors(9,:), '--'); 

[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(end), x(end), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(5), x(5), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop_02', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'HB');
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

[LI, ~, ~,x] = find_special_points('run_symm', 'HB');
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

% [LI, ~, ~,x] = find_special_points('run_zero', 'HB');
% scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
% scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(4), x(4), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(4), x(4), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop_03', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

[LI, ~, ~,x] = find_special_points('run_loop_03', 'SN');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));

xlim([-0.25, 0.15]);
ylim([0.12, 1.01]);
box on;
xlabel('\lambda_i');
ylabel('R_1')
title('N = 16, g = 2.5');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g2.5_N16.eps', 'Resolution', 300);

%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_zero', colors(2,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');
plot_curve('run_loop_02', colors(5,:), '--');
plot_curve('run_loop_03', colors(5,:), '--');
%plot_curve_po('run_anti_HB_00', colors(9,:), '--'); 
%plot_curve_po('run_symm_HB_00', colors(9,:), '--'); 
%plot_curve_po('run_loop_HB_00', colors(9,:), '--'); 

[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(end), x(end), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(5), x(5), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop_02', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'HB');
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

[LI, ~, ~,x] = find_special_points('run_symm', 'HB');
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

% [LI, ~, ~,x] = find_special_points('run_zero', 'HB');
% scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
% scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
[LI, ~, ~,x] = find_special_points('run_symm', 'BP');
scatter(LI(4), x(4), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(4), x(4), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop_03', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

[LI, ~, ~,x] = find_special_points('run_loop_03', 'SN');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));
scatter(LI(2), x(2), 50, 'filled', 'd', 'MarkerEdgeColor', colors(8,:), 'MarkerFaceColor',colors(8,:));

xlim([-0.26,  -0.23]);
ylim([-0.01, 0.18]);
box on;
xlabel('\lambda_i');
ylabel('R_1')
title('N = 16, g = 2.5');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g2.5_N16_zoomed.eps', 'Resolution', 300);

