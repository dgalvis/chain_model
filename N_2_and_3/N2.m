
%LE = .01;
omega = 2*pi;
c = 10;
g = 10;

LE_list = linspace(-2, 2, 1001);
Jtriv = zeros(size(LE_list));
Jsync = zeros(size(LE_list));

for i = 1:length(LE_list)

    LE = LE_list(i);
%% Trivial solution
% J(0,0)
J = zeros(4,4);

% Self
J(1,1) = LE-g; %dF1/dx1
J(2,2) = LE-g; %dF2/dx2

J(3,3) = LE-g;%dG1/dy1
J(4,4) = LE-g;%dG2/dy2

% Coupling
J(1,2) = g;
J(2,1) = g;

J(3,4) = g;
J(4,3) = g;

% xy
J(1,3) = -(omega + LE*c);
J(3,1) = omega + LE*c;

J(2,4) = -(omega+LE*c);
J(4,2) = omega+LE*c;

J

Jtriv(i) = max(real(eig(J)));

%% Synchronous solution

%dR1/dt = LE*R1*(1-R1^2) + g * (R2 cos(phi) - R1)
%dR2/dt = LE*R2*(1-R2^2) + g * (R1 cos(phi) - R2)
%dphi/dt = LE*c(R1^2 - R2^2) - g [R1/R2 + R2/R1] sin(phi2)

%J(1,1,0)

J = zeros(3,3);

J(1,1) = LE-3*LE-g;
J(2,2) = LE-3*LE-g;
J(3,3) = - 2*g*cos(0);

J(1,2) = g*cos(0);
J(2,1) = g*cos(0);

J(1,3) = - g*sin(0);
J(2,3) = -g*sin(0);

J(3,1) = 2*LE*c;
J(3,2) = - 2*LE*c;

Jsync(i) = max(real(eig(J)));

end

figure();hold all;
plot(LE_list, Jtriv, 'k');
plot(LE_list, Jsync, 'b');
scatter(0,0,'filled');