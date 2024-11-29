%%
clear;clc;
addpath('../../func_aux');

%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');
plot_curve('run_loop_02', colors(5,:), '--');
%plot_curve_po('run_anti_HB_00', colors(9,:), '--');                             

[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(end), x(end), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_loop_02', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

xlim([-1, 1]);
ylim([0.6, 1.01]);
box on;
xlabel('\lambda_i');
ylabel('R_1')
title('N = 10, g = 0.9');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g0.9_N10.eps', 'Resolution', 300);


%%
figure();hold all;
colors = get_colors();


plot_curve('run_sync', colors(1,:), '--');
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');
plot_curve('run_loop_02', colors(5,:), '--');
%plot_curve_po('run_anti_HB_00', colors(9,:), '--');                             

[LI, ~, ~,x] = find_special_points('run_sync', 'BP');
scatter(LI(end), x(end), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));


[LI, ~, ~,x] = find_special_points('run_loop', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_loop_02', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
scatter(LI(2), x(2), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'HB');
scatter(LI(1), x(1), 50, 'filled', 'd', 'MarkerEdgeColor', colors(7,:), 'MarkerFaceColor',colors(7,:));

xlim([-0.5505, -0.546]);
ylim([0.639, 0.641]);
box on;
xlabel('\lambda_i');
ylabel('R_1')
title('N = 10, g = 0.9, inset');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, '../g0.9_N10_zoom.eps', 'Resolution', 300);

