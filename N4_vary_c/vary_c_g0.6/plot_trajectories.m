%%
clear;clc;close all;
addpath('../../func_aux');
%%
N = 4;

bd = coco_bd_read('run_symm_c');
idx = find_idx(bd);
[LI, gg, cc, x, stab] = find_curves_all_x(bd, idx, N);

[~, IDX1] = min(abs(cc - 1));
[~, IDX0] = min(abs(cc - 0));
[~, IDXn1] = min(abs(cc - (-1)));


paux = [1, LI(IDX1), gg(IDX1), cc(IDX1), 2*pi];
plot_graph(x(:, IDX1)', paux, N);
exportgraphics(gcf, 'symm_c1.eps', 'Resolution', 300);

paux = [1, LI(IDX0), gg(IDX0), cc(IDX0), 2*pi];
plot_graph(x(:, IDX0)', paux, N);
exportgraphics(gcf, 'symm_c0.eps', 'Resolution', 300);

paux = [1, LI(IDXn1), gg(IDXn1), cc(IDXn1), 2*pi];
plot_graph(x(:, IDXn1)', paux, N);
exportgraphics(gcf, 'symm_cn1.eps', 'Resolution', 300);

%%
N = 4;

bd = coco_bd_read('run_anti_c');
idx = find_idx(bd);
[LI, gg, cc, x, stab] = find_curves_all_x(bd, idx, N);

[~, IDX1] = min(abs(cc - 1));
[~, IDX0] = min(abs(cc - 0));
[~, IDXn1] = min(abs(cc - (-1)));


paux = [1, LI(IDX1), gg(IDX1), cc(IDX1), 2*pi];
plot_graph(x(:, IDX1)', paux, N);
exportgraphics(gcf, 'anti_c1.eps', 'Resolution', 300);

paux = [1, LI(IDX0), gg(IDX0), cc(IDX0), 2*pi];
plot_graph(x(:, IDX0)', paux, N);
exportgraphics(gcf, 'anti_c0.eps', 'Resolution', 300);

paux = [1, LI(IDXn1), gg(IDXn1), cc(IDXn1), 2*pi];
plot_graph(x(:, IDXn1)', paux, N);
exportgraphics(gcf, 'anti_cn1.eps', 'Resolution', 300);


%%
addpath('../../functions');
N = 4;
LE = 1;
LI = -1;
g  = 0.6;
c  = 0;
omega = 2*pi;
p = [LE; LI; g; c; omega];

rng(1)
paux = p;
y0 = rand(1, 2*N-1);
tspan = [0,1000000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);


plot_graph(y(end, :), p, N);
exportgraphics(gcf, 'symm_c0.eps', 'Resolution', 300);
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
