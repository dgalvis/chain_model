%% Initialise
clear;clc;close all;
addpath('../functions');

%% N = 5

rng(1);
Nb = 5;
LEb = 1;
LIb = -0.4;
gb = 0.9;
cb  = 1;
omegab = 2*pi;
pb = [LEb; LIb; gb; cb; omegab];


y0 = rand(1, 2*Nb);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, pb, Nb));
[~,y] = ode15s(@(t,y)SL_vf(t, y, pb, Nb), tspan, y0, options);

y_final = y(end,:);
R = sqrt(y_final(1:Nb).^2 + y_final(Nb+1:end).^2);
T = atan2(y_final(Nb+1:end), y_final(1:Nb));
T = T - T(1);
T((Nb+1)/2) = nan; % only reasonable if R -> 0

R((Nb+1)/2) % needs to be ~0
plot_graph([R,T(2:end)], pb, Nb);pause(0.1);

%% N = 9

rng(1);
Nb = 9;
LEb = 1;
LIb = -0.1;
gb = 0.9;
cb  = 1;
omegab = 2*pi;
pb = [LEb; LIb; gb; cb; omegab];


y0 = rand(1, 2*Nb);
tspan = [0,1000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, pb, Nb));
[~,y] = ode15s(@(t,y)SL_vf(t, y, pb, Nb), tspan, y0, options);

y_final = y(end,:);
R = sqrt(y_final(1:Nb).^2 + y_final(Nb+1:end).^2);
T = atan2(y_final(Nb+1:end), y_final(1:Nb));
T = T - T(1);
T((Nb+1)/2) = nan;

R((Nb+1)/2)
plot_graph([R,T(2:end)], pb, Nb);pause(0.1);

abs(nanmean(R.*(cos(T) + 1j*sin(T))))
abs(nanmean((cos(T) + 1j*sin(T))))

%% N = 11

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

y_final = y(end,:);
R = sqrt(y_final(1:Nb).^2 + y_final(Nb+1:end).^2);
T = atan2(y_final(Nb+1:end), y_final(1:Nb));
T = T - T(1);
T((Nb+1)/2) = nan;

R((Nb+1)/2)
plot_graph([R,T(2:end)], pb, Nb);pause(0.1);

%% N = 21
rng(1);
Nb = 21;
LEb = 1;
LIb = -0.02;
gb = 0.9;
cb  = 1;
omegab = 2*pi;
pb = [LEb; LIb; gb; cb; omegab];


y0 = rand(1, 2*Nb);
tspan = [0,2000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, pb, Nb));
[~,y] = ode15s(@(t,y)SL_vf(t, y, pb, Nb), tspan, y0, options);

y_final = y(end,:);
R = sqrt(y_final(1:Nb).^2 + y_final(Nb+1:end).^2);
T = atan2(y_final(Nb+1:end), y_final(1:Nb));
T = T - T(1);
T((Nb+1)/2) = nan;


R((Nb+1)/2)
plot_graph([R,T(2:end)], pb, Nb);

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

