function [J, DFDP] = SL_polar_jac(~, y, p, N, varargin) 

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
    %omega   = p(5,:);
    

    R = y(1:N, :);
    theta = zeros(size(R));

    if strcmp(flag, 'moving_frame') 
        theta(2:N, :) = wrapToPi(y(N+1:end, :));
    else
        theta = wrapToPi(y(N+1:end, :));
    end


    J = zeros(2*size(R, 1), 2*size(R,1), size(R,2));
    JRR = zeros(size(R, 1), size(R,1), size(R,2));
    JRT = zeros(size(R, 1), size(R,1), size(R,2));
    JTR = zeros(size(R, 1), size(R,1), size(R,2));
    JTT = zeros(size(R, 1), size(R,1), size(R,2));

    % JRR
    aux1 = zeros(size(R));
    aux1(1,:) = lambdaE - g - 3*lambdaE.*R(1,:).^2;
    aux1(2:N-1,:) = lambdaI - 2*g - 3*lambdaI.*R(2:N-1,:).^2;
    aux1(N, :) = lambdaE - g - 3*lambdaE.*R(N,:).^2;

    aux2 = g.*cos(theta(2:N,:) - theta(1:N-1,:));% R upper off diagonal
    aux3 = g.*cos(theta(1:N-1,:) - theta(2:N,:));% R lower off diagonal

    for i = 1:N
        JRR(i, i, :) = aux1(i, :);
        if i~= N
            JRR(i, i+1, :) = aux2(i, :);
        end
        if i~=1
            JRR(i, i-1, :) = aux3(i-1, :);
        end
    end


    % JRT
    aux1 = zeros(size(R));
    aux1(1,:)     = g.*R(2,:).*sin(theta(2,:) - theta(1,:));
    aux1(2:N-1,:) = g.*R(3:N,:).*sin(theta(3:N,:) - theta(2:N-1,:)) + g.*R(1:N-2,:).*sin(theta(1:N-2,:) - theta(2:N-1,:));
    aux1(N,:)     = g.*R(N-1,:).*sin(theta(N-1,:) - theta(N,:));

    aux2 = -g.*R(2:N,:).*sin(theta(2:N,:) - theta(1:N-1,:));% 
    aux3 = -g.*R(1:N-1,:).*sin(theta(1:N-1,:) - theta(2:N,:));% 

    for i = 1:N
        JRT(i, i, :) = aux1(i, :);
        if i~= N
            JRT(i, i+1, :) = aux2(i, :);
        end
        if i~=1
            JRT(i, i-1, :) = aux3(i-1, :);
        end
    end    

    % JTR
    aux1 = zeros(size(R));
    aux1(1,:)     = -2.*c.*lambdaE.*R(1,:)     - g.*(R(2,:)./R(1,:).^2.*sin(theta(2,:) - theta(1,:)));
    aux1(2:N-1,:) = -2.*c.*lambdaI.*R(2:N-1,:) - g.*(R(1:N-2,:)./R(2:N-1,:).^2.*sin(theta(1:N-2,:) - theta(2:N-1,:)) + R(3:N,:)./R(2:N-1,:).^2.*sin(theta(3:N,:) - theta(2:N-1,:)));
    aux1(N,:)     = -2.*c.*lambdaE.*R(N,:)     - g.*(R(N-1,:)./R(N,:).^2.*sin(theta(N-1,:) - theta(N,:)));


    aux2 = g./R(1:N-1,:).*sin(theta(2:N,:) - theta(1:N-1,:));% 
    aux3 = g./R(2:N,:).*sin(theta(1:N-1,:) - theta(2:N,:));% 

    for i = 1:N
        JTR(i, i, :) = aux1(i, :);
        if i~= N
            JTR(i, i+1, :) = aux2(i, :);
        end
        if i~=1
            JTR(i, i-1, :) = aux3(i-1, :);
        end
    end   

    % JTT
    aux1 = zeros(size(R));
    aux1(1,:) = -g.*(R(2,:)./R(1,:).*cos(theta(2,:) - theta(1,:)));
    aux1(2:N-1,:) = -g.*(R(1:N-2,:)./R(2:N-1,:).*cos(theta(1:N-2,:) - theta(2:N-1,:)) + R(3:N,:)./R(2:N-1,:).*cos(theta(3:N,:) - theta(2:N-1,:)));
    aux1(N,:) = -g.*(R(N-1,:)./R(N,:).*cos(theta(N-1,:) - theta(N,:)));  

    aux2 = g.*R(2:N,:)./R(1:N-1,:).*cos(theta(2:N,:) - theta(1:N-1,:));% 
    aux3 = g.*R(1:N-1,:)./R(2:N,:).*cos(theta(1:N-1,:) - theta(2:N,:));% 

    for i = 1:N
        JTT(i, i, :) = aux1(i, :);
        if i~= N
            JTT(i, i+1, :) = aux2(i, :);
        end
        if i~=1
            JTT(i, i-1, :) = aux3(i-1, :);
        end
    end   


    if strcmp(flag, 'standard_frame')
        J(1:N, 1:N, :) = JRR;
        J(1:N, N+1:end, :) = JRT;
        J(N+1:end, 1:N,:) = JTR;
        J(N+1:end, N+1:end,:) = JTT;
    elseif strcmp(flag, 'moving_frame')
        J(1:N, 1:N, :) = JRR;
        J(1:N, N+1:end, :) = JRT;
        J(N+1:end, 1:N,:) = JTR;
        J(N+1:end, N+1:end,:) = JTT;
        J(N+1:end, :,:) = J(N+1:end, :,:) - J(N+1,:,:);  
        J = J([1:N, N+2:end], [1:N, N+2:end],:);
    else
        error('invalid flag');
    end
   
    if nargout == 2
        DFDP = SL_polar_DFDP(0, y, p, N, 'flag', flag);
    end
end



% % Testing Code to compare moving and standard frame
% clear;clc;close all;
% addpath('functions');
% rng('shuffle');
% N = 10;
% %lambdaE, lambdaI, g, c, omega
% p = [1, -1, 0.01*(N-1)^2, 1, 2*pi]';
% 
% 
% options_ode = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Vectorized', 'on', 'Jacobian', @(t,y)SL_polar_jac(t,y,p,N, 'flag', 'standard_frame'));
% y0    = rand(2*N, 1);
% tspan = 0:0.01:100;
% 
% [t_1, y_1] = ode15s(@(t,y)SL_polar_vf(t, y, p, N, 'flag', 'standard_frame'), tspan, y0, options_ode);
% 
% figure();hold all;
% plot(t_1, sum(y_1(:, 1:N).^2,2)/N);
% xlabel('time');
% ylabel('L^2/N');
% 
% 
% 
% options_ode = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Vectorized', 'on', 'Jacobian', @(t,y)SL_polar_jac(t,y,p,N, 'flag', 'moving_frame'));
% 
% y0(N+1:end) = y0(N+1:end) - y0(N+1);
% y0(N+1) = [];
% [t_2, y_2] = ode15s(@(t,y)SL_polar_vf(t, y, p, N, 'flag', 'moving_frame'), tspan, y0, options_ode);
% 
% y_1(:, N+1:end) = y_1(:, N+1:end) - y_1(:, N+1);
% y_1(:, N+1) = [];
% max(abs(y_2 - y_1))