% The vector field is vectorized
% Moving frame phi_i = theta_i - theta_1 (2*N-1 variables since theta_1 = 0)
% Standard frame (2*N variables)
% Zeros only (in moving frame) can aid in numerical continuation of the fixed points because
% it avoids the problem of divide by zero. This cannot be used for
% simulation or for identifying the stability of those fixed points
function [dydt, dtheta] = SL_polar_vf(~, y, p, N, varargin)

    flag = 'moving_frame';
    for i = 1:length(varargin)
        if strcmp(varargin{i}, 'flag')
            flag = varargin{i+1};
        end
    end


    lambdaE = p(1,:);
    lambdaI = p(2,:);
    g       = p(3,:);
    c       = p(4,:);
    omega   = p(5,:);
    

    R = y(1:N, :);
    theta = zeros(size(R));

    if strcmp(flag, 'moving_frame') || strcmp(flag, 'zeros_only')
        theta(2:N, :) = wrapToPi(y(N+1:end, :));
    else
        theta = wrapToPi(y(N+1:end, :));
    end

    dRdt = zeros(N, size(y,2));
    dtheta = zeros(N, size(y,2));

    dRdt(1,:)        = lambdaE.*R(1, :).*(1-R(1, :).^2);
    dRdt(2:N-1, :)   = lambdaI.*R(2:N-1, :).*(1-R(2:N-1, :).^2);
    dRdt(N, :)       = lambdaE.*R(N, :).*(1-R(N, :).^2);

    dtheta(1, :)     = omega + c.*lambdaE.*(1-R(1, :).^2);
    dtheta(2:N-1, :) = omega + c.*lambdaI.*(1-R(2:N-1,:).^2);
    dtheta(N, :)     = omega + c.*lambdaE.*(1-R(N, :).^2);


    dRdt(1,:)        = dRdt(1,:)      + g.*(R(2, :).*cos(theta(2, :) - theta(1, :)) - R(1, :));
    dRdt(2:N-1, :)   = dRdt(2:N-1, :) + g.*(R(1:N-2, :).*cos(theta(1:N-2, :) - theta(2:N-1, :)) + R(3:N, :).*cos(theta(3:N, :) - theta(2:N-1, :)) - 2*R(2:N-1, :));
    dRdt(N, :)       = dRdt(N, :)     + g.*(R(N-1, :).*cos(theta(N-1, :) - theta(N, :)) - R(N, :));


    if strcmp(flag, 'moving_frame') || strcmp(flag, 'standard_frame')
        % For simulation
        dtheta(1, :)     = dtheta(1, :)     + g.*(R(2, :)./R(1, :).*sin(theta(2, :) - theta(1, :)));
        dtheta(2:N-1, :) = dtheta(2:N-1, :) + g.*(R(1:N-2, :)./R(2:N-1, :).*sin(theta(1:N-2, :) - theta(2:N-1, :)) + R(3:N, :)./R(2:N-1, :).*sin(theta(3:N, :) - theta(2:N-1, :)));
        dtheta(N, :)     = dtheta(N, :)     + g.*(R(N-1, :)./R(N, :).*sin(theta(N-1, :) - theta(N, :)));

        if strcmp(flag, 'moving_frame')
            dydt = [dRdt; dtheta(2:end,:) - dtheta(1, :)];
        elseif strcmp(flag, 'standard_frame')
            dydt = [dRdt; dtheta];
        end
    elseif strcmp(flag, 'zeros_only')
        % Not for simulation, only for finding zeros!!!! (gets rid of
        % divide by zero issues)
        dtheta(1, :)     = dtheta(1, :).*R(1,:)     + g.*(R(2, :).*sin(theta(2, :) - theta(1, :)));
        dtheta(2:N-1, :) = dtheta(2:N-1, :).*R(2:N-1,:) + g.*(R(1:N-2, :).*sin(theta(1:N-2, :) - theta(2:N-1, :)) + R(3:N, :).*sin(theta(3:N, :) - theta(2:N-1, :)));
        dtheta(N, :)     = dtheta(N, :).*R(N,:)     + g.*(R(N-1, :).*sin(theta(N-1, :) - theta(N, :)));
        dydt = [dRdt; dtheta(2:end,:).*R(1,:) - dtheta(1, :).*R(2:N,:)];
    else
        error('flag is invalid')
    end
end

% Testing Code to compare moving and standard frame
%clear;clc;close all;
%addpath('functions');
%rng('shuffle');
%N = 10;
% % lambdaE, lambdaI, g, c, omega
%p = [1, -1, 0.01*(N-1)^2, 1, 2*pi]';


%options_ode = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Vectorized', 'on');%, 'Jacobian', @(t,y)SL_polar_jac(t,y,p,N, true));

%y0    = rand(2*N, 1);
%tspan = 0:0.01:100;

%[t_1, y_1] = ode15s(@(t,y)SL_polar_vf(t, y, p, N, 'flag', 'standard_frame'), tspan, y0, options_ode);

% figure();hold all;
% 
% plot(t_1, sum(y_1(:, 1:N).^2,2)/N);
% 
% xlabel('time');
% ylabel('L^2/N');
% 
% y0(N+1:end) = y0(N+1:end) - y0(N+1);
% y0(N+1) = [];
% [t_2, y_2] = ode15s(@(t,y)SL_polar_vf(t, y, p, N, 'flag', 'moving_frame'), tspan, y0, options_ode);
% 
% y_1(:, N+1:end) = y_1(:, N+1:end) - y_1(:, N+1);
% y_1(:, N+1) = [];
% max(abs(y_2 - y_1))
