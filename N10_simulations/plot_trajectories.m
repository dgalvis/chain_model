%% Initialise
clear;clc;close all;
addpath('../functions');
%% Baseline Parameters

N = 10;
LE = 1;
LI = 1;
g  = 0.9;
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];


%% Sync
rng(1);
paux = p;
paux(2) = 1;

y0 = [ones(1,N), zeros(1, N-1)];
tspan = [0,10000000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

max(abs(SL_polar_vf(0, y(end,:)', paux, N)))

plot_graph(y, paux, N);

exportgraphics(gcf, 'sync.eps', 'Resolution', 300);

%% Symmetric Wave
rng(13);
paux = p;
paux(2) = -0.06;

y0 = rand(1, 2*N-1);
tspan = [0,1000000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

max(abs(SL_polar_vf(0, y(end,:)', paux, N)))

plot_graph(y, paux, N);

exportgraphics(gcf, 'symm.eps', 'Resolution', 300);

%% Anti-Symmetric Wave
rng(1);
paux = p;
paux(2) = -0.06;

y0 = rand(1, 2*N-1);
tspan = [0,1000000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

max(abs(SL_polar_vf(0, y(end,:)', paux, N)))

plot_graph(y, paux, N);
exportgraphics(gcf, 'anti.eps', 'Resolution', 300);


%% EXTRA N = 11

rng(1);
Nb = 11;
LEb = 1;
LIb = -0.06;
gb = 0.9;
cb  = 1;
omegab = 2*pi;
pb = [LEb; LIb; gb; cb; omegab];


y0 = rand(1, 2*Nb);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, pb, Nb));
[~,y] = ode15s(@(t,y)SL_vf(t, y, pb, Nb), tspan, y0, options);

%%
y_final = y(end,:);

R = sqrt(y_final(1:Nb).^2 + y_final(Nb+1:end).^2);
T = atan2(y_final(Nb+1:end), y_final(1:Nb));
T = T - T(1);
T((Nb+1)/2) = nan;

plot_graph([R,T(2:end)], pb, Nb);
exportgraphics(gcf, 'anti_N11.eps', 'Resolution', 300);

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
%shading flat
set(gca, 'YDir', 'reverse');
ylabel('T (periods)');

%title('theta_i vs. t');
cobj = colorbar;
cobj.Ticks = [-pi, 0, pi];
cobj.TickLabels = {'-\pi', '0', '\pi'};
clim([-pi, pi]);

xticks(1:N);
yticks(linspace(tidx(1), tidx(end), 6));
yticklabels(0:5);
set(gca, 'Fontsize', 20);

nexttile;hold all;
plot(1:N, R, 'k');
scatter(1:N, R, 100, 'filled', 'k');
ylim([0,1]);
xlim([0.5, N+0.5]);
set(gca, 'YDir', 'normal');
ylabel('R_j');
xlabel('N');
xticks(1:N);

%sgtitle(['N=', num2str(N), 'g=', num2str(p(3)), 'LI=', num2str(p(2))]);
set(gca, 'Fontsize', 20);
set(gcf, 'Position', [100, 100, 400, 400]);

end

