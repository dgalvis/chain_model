%%
clear;clc;close all;
addpath('../../func_aux')
%%
figure();hold all;
colors = get_colors();
plot_curve('run_symm', colors(3,:), '--');
plot_curve('run_anti', colors(4,:), '--');
plot_curve('run_loop', colors(5,:), '--');

[LI, ~, ~, x] = find_special_points('run_symm', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

[LI, ~, ~,x] = find_special_points('run_anti', 'BP');
scatter(LI(1), x(1), 50, 'filled', 's', 'MarkerEdgeColor', colors(6,:), 'MarkerFaceColor',colors(6,:));

xlim([-0.1, 0.02]);
ylim([0.96, 0.977]);
xlabel('\lambda_i');
ylabel('R_1');
box on;
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties


%%
N = 4;
IDX0 = 1;
IDX1 = 201;
IDX2 = 401;
IDX3 = 601;
IDX4 = 801;
IDX5 = 490;

%%


bd = coco_bd_read('run_loop');
idx = find_idx(bd);
[LI, ~, ~, x, ~] = find_curves_all_x(bd, idx, N);

scatter(LI(IDX0), x(1, IDX0), 100, 'filled', 'bs');
scatter(LI(IDX1), x(1, IDX1), 100, 'filled', 'k^')
scatter(LI(IDX2), x(1, IDX2), 100, 'filled', 'ko')
scatter(LI(IDX3), x(1, IDX3), 100, 'filled', 'ks')
scatter(LI(IDX4), x(1, IDX4), 100, 'filled', 'kv')
scatter(LI(IDX5), x(1, IDX5), 100, 'filled', 'bo');

exportgraphics(gcf, 'loop_R1.eps', 'Resolution', 300);

%%
figure();hold all;
plot(LI, x(end, :),  '--', 'color', colors(5,:));


% Add in sym/anti curves
bds = coco_bd_read('run_symm');
bda = coco_bd_read('run_anti');

idxs = find_idx(bds);
idxa = find_idx(bda);

[LIs, ~, ~, xs, stabs] = find_curves_all_x(bds, idxs, N);
[LIa, ~, ~, xa, staba] = find_curves_all_x(bda, idxa, N);

[stabs, ustabs] = find_intervals(stabs);
[staba, ustaba] = find_intervals(staba);

plot(LIs(ustabs(1,1):ustabs(1,2)), xs(end, ustabs(1,1):ustabs(1,2)),  '--', 'color', colors(3,:));
plot(LIa(ustaba(1,1):ustaba(1,2)), xa(end, ustaba(1,1):ustaba(1,2)),  '--', 'color', colors(4,:));
plot(LIs(stabs(1,1):stabs(1,2)), xs(end, stabs(1,1):stabs(1,2)),  '-', 'color', colors(3,:));
plot(LIa(staba(1,1):staba(1,2)), xa(end, staba(1,1):staba(1,2)),  '-', 'color', colors(4,:));

% add in phi_N + 2pi
plot(LIs(ustabs(1,1):ustabs(1,2)), xs(end, ustabs(1,1):ustabs(1,2)) + 2*pi,  '--', 'color', colors(3,:));
plot(LIs(stabs(1,1):stabs(1,2)), xs(end, stabs(1,1):stabs(1,2)) + 2*pi,  '-', 'color', colors(3,:));


axis tight;
xlim([-0.11, 0.02]);
ylim([0-0.1, 2*pi+0.1]);
xlabel('\lambda_i');
ylabel('\phi_N');
box on;
set(gca, 'FontSize', 16, 'LineWidth', 1); % Customize axis properties

%%
N = 4;

bd = coco_bd_read('run_loop');
idx = find_idx(bd);
[LI, gg, cc, x, stab] = find_curves_all_x(bd, idx, N);

scatter(LI(IDX0), x(end, IDX0), 100, 'filled', 'bs');
scatter(LI(IDX0), x(end, IDX0) + 2*pi, 100, 'filled', 'bs');

scatter(LI(IDX1), x(end, IDX1), 100, 'filled', 'k^');
scatter(LI(IDX2), x(end, IDX2), 100, 'filled', 'ko');
scatter(LI(IDX3), x(end, IDX3), 100, 'filled', 'ks');
scatter(LI(IDX4), x(end, IDX4), 100, 'filled', 'kv');
scatter(LI(IDX5), x(end, IDX5), 100, 'filled', 'bo');

exportgraphics(gcf, 'loop_phiN.eps', 'Resolution', 300);

%%
LE = 1;
omega = 2*pi;
g = 0.1;
c = 1.0;

p = [LE, LI(IDX0), g, c, omega];
plot_graph(x(:, IDX0)', p, N)
exportgraphics(gcf, 'blue_square.eps', 'Resolution', 300);

p = [LE, LI(IDX1), g, c, omega];
plot_graph(x(:, IDX1)', p, N)
exportgraphics(gcf, 'up.eps', 'Resolution', 300);


p = [LE, LI(IDX2), g, c, omega];
plot_graph(x(:, IDX2)', p, N)
exportgraphics(gcf, 'circle.eps', 'Resolution', 300);


p = [LE, LI(IDX3), g, c, omega];
plot_graph(x(:, IDX3)', p, N)
exportgraphics(gcf, 'square.eps', 'Resolution', 300);


p = [LE, LI(IDX4), g, c, omega];
plot_graph(x(:, IDX4)', p, N);
exportgraphics(gcf, 'down.eps', 'Resolution', 300);

p = [LE, LI(IDX5), g, c, omega];
plot_graph(x(:, IDX5)', p, N)
exportgraphics(gcf, 'blue_circle.eps', 'Resolution', 300);

%% Functions
function plot_graph(y, p, N)

LE = p(1);
LI = p(2);
g  = p(3);
c  = p(4);
omega = p(5);

R = y(end, 1:N);
T = [0, y(end, N+1:end)];
om = omega + c.*LE.*(1-R(1).^2) + g.*(R(2)./R(1).*sin(T(2) - T(1)));

sc = 2*pi/om;
t = (0:0.01:5*sc)';
tidx = 1:length(t);
A = R.*exp(1j*(om .* t + T)) ;


figure();tiledlayout(4, 1);
nexttile([3,1]);hold all;
imagesc(angle(A));
axis tight;
set(gca, 'YDir', 'reverse');
ylabel('T (periods)');

title('theta_i vs. t');
cobj = colorbar;
cobj.Ticks = [-pi, 0, pi];
cobj.TickLabels = {'-\pi', '0', '\pi'};
clim([-pi, pi]);

xticks(1:N);
yticks(linspace(tidx(1), tidx(end), 6));
yticklabels(0:5);

nexttile;hold all;
plot(1:N, R, 'k');
scatter(1:N, R, 100, 'filled', 'k');
ylim([0,1]);
xlim([0.5, N+0.5]);
set(gca, 'YDir', 'normal');
ylabel('R_i');
xlabel('N');
sgtitle(['N=', num2str(N), 'g=', num2str(p(3)), 'LI=', num2str(p(2))]);
xticks(1:N)

set(gcf, 'Position', [100, 100, 200, 400]);

end