%%
clear;clc;close all;
addpath('../../func_aux');
colors = get_colors();
glist = 1:0.1:2.5;
%figure();hold all;

gs_all = [];
xs_all = [];
LIs_all = [];
ga_all = [];
xa_all = [];
LIa_all = [];
gh_all = [];
xh_all = [];
LIh_all = [];
for i = 1:length(glist)

    g = glist(i);
    [LIs, ~, ~,xs] = find_special_points(['run_symm_',num2str(i)], 'BP');
    [LIa, ~, ~,xa] = find_special_points(['run_anti_',num2str(i)], 'BP');

    [LIh, ~, ~,xh] = find_special_points(['run_anti_',num2str(i)], 'HB');

    LIs = LIs(xs > 0.1& xs < 0.995);
    xs = xs(xs > 0.1& xs < 0.995);

    LIa = LIa(xa > 0.1 & xa < 0.995);
    xa = xa(xa > 0.1& xa < 0.995);

    LIh = LIh(xh > 0.1 & xh < 0.995);
    xh = xh(xh > 0.1& xh < 0.995);

    gs = g*ones(size(LIs));
    ga = g*ones(size(LIa));
    gh = g*ones(size(LIh));

    gs_all = [gs_all; gs];
    LIs_all = [LIs_all; LIs];
    xs_all = [xs_all; xs];

    ga_all = [ga_all; ga];
    LIa_all = [LIa_all; LIa];
    xa_all = [xa_all; xa];

    gh_all = [gh_all; gh];
    LIh_all = [LIh_all; LIh];
    xh_all = [xh_all; xh];

   %scatter3(g*ones(size(LIs)), LIs, xs, 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
   %scatter3(g*ones(size(LIa)), LIa, xa, 50, 'filled', 'o', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));
end


%%

figure(); hold on;



% Shade the area between the first set of curves
fill([LIs_all(:, 1); flipud(LIa_all(:, 1))], ...
     [gs_all(:, 1); flipud(ga_all(:, 1))], ...
     [0.8, 0.8, 0.8], 'EdgeColor', 'none'); % Light red shading

% Plot the first set of curves
plot(LIs_all(:, 1), gs_all(:, 1), 'b', 'LineWidth', 0.5);
plot(LIa_all(:, 1), ga_all(:, 1), 'r', 'LineWidth', 0.5);


% Shade the area between the second set of curves
fill([LIs_all(:, 2); flipud(LIa_all(:, 2))], ...
     [gs_all(:, 2); flipud(ga_all(:, 2))], ...
     [0.7, 0.7, 0.7], 'EdgeColor', 'none'); % Slightly darker red shading


% Plot the second set of curves
plot(LIs_all(:, 2), gs_all(:, 2), 'b', 'LineWidth', 1);
plot(LIa_all(:, 2), ga_all(:, 2), 'r', 'LineWidth', 1);


% Shade the area between the third set of curves
fill([LIs_all(:, 3); flipud(LIh_all(:, 1))], ...
     [gs_all(:, 3); flipud(gh_all(:, 1))], ...
     [0.4, 0.4, 0.4], 'EdgeColor', 'none'); % Darker red shading

% Shade the area between the third set of curves
fill([LIh_all(:, 1); flipud(LIa_all(:, 3))], ...
     [gh_all(:, 1); flipud(ga_all(:, 3))], ...
     colors(3,:), 'EdgeColor', 'none'); % Darker red shading

% Plot the third set of curves
plot(LIs_all(:, 3), gs_all(:, 3), 'b', 'LineWidth', 1);
plot(LIa_all(:, 3), ga_all(:, 3), 'r', 'LineWidth', 1);

plot(LIh_all, gh_all, 'k--', 'linewidth', 2);

% Customize the axes and appearance
axis tight;
box on;
xlabel('\lambda_i', 'Interpreter', 'tex');
ylabel('g', 'Interpreter', 'tex');
zlabel('R_1', 'Interpreter', 'tex');
title('N = 26');
set(gca, 'FontSize', 16, 'LineWidth', 1);

% Save the figure
exportgraphics(gcf, '../N26_2par.eps', 'Resolution', 300);
