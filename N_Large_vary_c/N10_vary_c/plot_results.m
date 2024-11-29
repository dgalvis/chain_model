%%
clear;clc;
addpath('../../func_aux');
%% Original Image
% figure();hold all;
% thm = struct();
% thm.special = {'BP'};
% %thm.lspec = {'g-', 'LineWidth', 1};
% coco_plot_bd(thm, 'run_sync_c', 'lambdaI', 'x');
% 
% figure();hold all;
% coco_plot_bd(thm, 'run_symm_c', 'lambdaI', 'x');
% figure();hold all;
% coco_plot_bd(thm, 'run_anti_c', 'lambdaI', 'x');
% 
% box on;
% 
% title('g = 0.9  / N = 10');
% ylabel('R_1^*');


%%
figure();hold all;
colors = get_colors();


%plot_curve('run_sync_c', colors(1,:), '--', 'cc');
plot_curve('run_symm_c', colors(3,:), '--', 'cc');
plot_curve('run_anti_c', colors(4,:), '--', 'cc');                         
plot_curve('run_loop_c_1', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_2', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_3', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_4', colors(5,:), '--', 'cc');    


[~, ~, cc, x] = find_special_points('run_loop_c_1', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_2', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_3', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_4', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

xlim([-1.5, 1.5]);
box on;
xlabel('c');
ylabel('R_1')
title('g = 0.9, \lambda_i = -0.2');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, 'g0.9_N10_LIn0.2_vary_c.eps', 'Resolution', 300);

%%
figure();hold all;
colors = get_colors();


%plot_curve('run_sync_c', colors(1,:), '--', 'cc');
plot_curve('run_symm_c', colors(3,:), '--', 'cc');
plot_curve('run_anti_c', colors(4,:), '--', 'cc');                         
plot_curve('run_loop_c_1', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_2', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_3', colors(5,:), '--', 'cc');    
plot_curve('run_loop_c_4', colors(5,:), '--', 'cc');    


[~, ~, cc, x] = find_special_points('run_loop_c_1', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_2', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_3', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

[~, ~, cc, x] = find_special_points('run_loop_c_4', 'BP');
for i = 1:length(cc)
    scatter(cc(i), x(i), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end

xlim([1.255, 1.285]);
ylim([0.6995, 0.7050])
box on;
xlabel('c');
ylabel('R_1')
title('g = 0.9, \lambda_i = -0.2');
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties
exportgraphics(gcf, 'g0.9_N10_LIn0.2_vary_c_zoomed.eps', 'Resolution', 300);