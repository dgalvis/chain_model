%% Initialise
clear;clc;close all;
addpath('../functions');
%% Baseline Parameters

N = 4;
LE = 1;
LI = 1;
g  = 1;
c  = 1;
omega = 2*pi;
p = [LE; LI; g; c; omega];


%% Sync
rng(1)
paux = p;
paux(2) = 1;
y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

plot_graph(y, paux, N);

exportgraphics(gcf, 'sync.eps', 'Resolution', 300);

%% Symmetric Wave

rng(1);
paux = p;
paux(2) = -1;

y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

plot_graph(y, paux, N);

exportgraphics(gcf, 'symm.eps', 'Resolution', 300);


%% Anti-Symmetric Wave
rng(7);
paux = p;
paux(2) = -1;

y0 = rand(1, 2*N-1);
tspan = [0,10000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

plot_graph(y, paux, N);

exportgraphics(gcf, 'anti.eps', 'Resolution', 300);


%% Asymmetric Loop Stable
rng(7);
paux = p;
paux(2) = -2.2;
paux(3) = 1.12;

y0 = rand(1, 2*N-1);
tspan = [0,100000];
options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
[~,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);

plot_graph(y, paux, N);

exportgraphics(gcf, 'asym.eps', 'Resolution', 300);

%% Torus solution
% rng(7);
% paux = p;
% paux(2) = -2;
% paux(3) = 1.12;
% 
% y0 = rand(1, 2*N-1);
% tspan = [0,2000];
% options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
% [t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);
% y0 = y(end,:);
% tspan = 0:0.01:800;
% [t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);
% 
% R = y(:,1:N);
% T = [zeros(length(y), 1), y(:, N+1:end)];
% 
% figure();
% plot(t,R, 'linewidth',2);
% legend('R_1', 'R_2', 'R_3', 'R_4');
% xlabel('time');
% ylabel('R (moving frame)');
% figure();
% plot(t,T, 'linewidth',2);
% legend('\phi_1', '\phi_2', '\phi_3', '\phi_4');
% xlabel('time');
% ylabel('phi (moving frame)');
% 
% 
% %% Torus solution part 2
% rng(1);
% paux = p;
% paux(2) = -2;
% paux(3) = 1.12;
% 
% y0 = rand(1, 2*N-1);
% tspan = [0,2000];
% options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_polar_jac(t,y, paux, N));
% [t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);
% y0 = y(end,:);
% tspan = 0:0.01:800;
% [t,y] = ode15s(@(t,y)SL_polar_vf(t, y, paux, N), tspan, y0, options);
% 
% R = y(:,1:N);
% T = [zeros(length(y), 1), y(:, N+1:end)];
% 
% figure();
% plot(t,R, 'linewidth',2);
% legend('R_1', 'R_2', 'R_3', 'R_4');
% xlabel('time');
% ylabel('R (moving frame)');
% figure();
% plot(t,T, 'linewidth',2);
% legend('\phi_1', '\phi_2', '\phi_3', '\phi_4');
% xlabel('time');
% ylabel('phi (moving frame)');
% 
% %% Simulate Torus further
% 
% y0 = y(end,:);
% R = y0(1:N);
% T = [0,y0(N+1:end)];
% 
% y0 = [R.*cos(T), R.*sin(T)];
% 
% options = odeset(AbsTol=1e-12, RelTol=1e-12, Jacobian=@(t,y)SL_jac(t,y, paux, N));
% tspan = 0:0.01:200;
% [t,y] = ode15s(@(t,y)SL_vf(t, y, paux, N), tspan, y0, options);
% 
% %% Torus again
% R = sqrt(y(:, 1:N).^2 + y(:, N+1:end).^2);
% T = atan2(y(:,N+1:end), y(:, 1:N));
% 
% y2 = [R, T];
% dxdt = SL_polar_vf(0, y2', paux, N, 'flag', 'standard_frame')';
% dTdt = dxdt(:, N+1:end);
% 
% figure();
% plot(t,R, 'linewidth',2);
% legend('R_1', 'R_2', 'R_3', 'R_4');
% xlabel('time');
% ylabel('R (moving frame)');
% 
% figure();
% plot(t,T, 'linewidth',2);
% legend('\theta_1', '\theta_2', '\theta_3', '\theta_4');
% xlabel('time');
% ylabel('theta (standard frame)');
% 
% figure();
% plot(t,dTdt, 'linewidth',2);
% legend('what_1', 'what_2', 'what_3', 'what_4');
% xlabel('time');
% ylabel('emergent frequency');


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

