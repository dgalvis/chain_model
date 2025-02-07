%% Initialise
clear; clc; close all;
addpath('../../functions'); % Add the 'functions' directory to the search path

%% Initial parameter values and Fsolve information

% Fixed parameters
LE = 1;                % Length of the element
c  = 2;                % Constant parameter
omega = 2*pi;          % Angular frequency

% Varied parameters
LI = nan;              % Initial value for the length parameter (will be varied)
g  = nan;              % Initial value for the gain parameter (will be varied)

% Parameter vector
p = [LE; LI; g; c; omega];

% Range of varied parameters
g_list = linspace(1e-6, 5, 2001);  % List of values for the gain parameter
N_list = [4, 10];                 % List of values for the number of elements

% Fsolve options
options = optimoptions('fsolve', 'MaxIterations', 10000, 'StepTolerance', 1e-12, 'Display', 'none');

%% Making sure the curve is monotonic before looking for the zero
% p(3) = 3;
% LI_list = linspace(-3, 0, 1001);
% N = 10;
% for i = 1:length(LI_list)
%     x(i) = trivial_zero_stab(LI_list(i), p, N);
% end
% plot(LI_list, x)


%% Main computation loop
triv_0_list = nan(length(N_list), length(g_list)); % Preallocate matrix for results of trivial_zero_stab
triv_1_list = nan(length(N_list), length(g_list)); % Preallocate matrix for results of trivial_one_stab

% Loop over the list of N values
for i = 1:length(N_list)
    % Loop over the list of g values
    for j = 1:length(g_list)
        % display iteration
        [i, j]

        % Set changing parameters
        N    = N_list(i);         % Current number of elements
        p(3) = g_list(j);         % Update parameter vector with current g value
      

        % CHECK MONOTONICITY (This is slow so comment out if not using
        % LI_list = linspace(-10, 0.1, 1001);
        % x = zeros(size(LI_list));
        % y = x;
        % for k = 1:length(LI_list)
        %     x(k) = trivial_zero_stab(LI_list(k), p, N);
        %     y(k) = trivial_one_stab(LI_list(k), p, N);
        % end
        % 
        % % Find zero crossings
        % sign_changes_x = diff(sign(x)); % Detect sign changes in x
        % sign_changes_y = diff(sign(y)); % Detect sign changes in y
        % 
        % zero_crossings_x = find(sign_changes_x ~= 0); % Indices of zero crossings 
        % zero_crossings_y = find(sign_changes_y ~= 0); % Indices of zero crossings
        % 
        % % Count the number of zero crossings
        % num_crossings_x = length(zero_crossings_x);
        % num_crossings_y = length(zero_crossings_y);
        % 
        % % The second statment is just checking that there isn't some degeneracy of the form:
        % % [pos, pos, 0, 0, 0] (as happens for g  = 0 with the synchronous solution)
        % if (num_crossings_x > 1) || (sum(sign(x) == 0) > 0) 
        %    [N, p(3), 'check x']
        %    pause
        % end
        % if (num_crossings_y > 1) || (sum(sign(y) == 0) > 0)
        %    [N, p(3), 'check y']
        %    pause;
        % end        
        % END CHECK MONOTONICITY

        % Solve for trivial_zero_stab
        [x0, fval0, exitflag0] = fsolve(@(x)trivial_zero_stab(x, p, N), -1, options);
        
        % Solve for trivial_one_stab
        [x1, fval1, exitflag1] = fsolve(@(x)trivial_one_stab(x, p, N), -1, options);
        
        % Store results if the solution converged
        if exitflag0 == 1
            triv_0_list(i, j) = x0;
        end
        if exitflag1 == 1
            triv_1_list(i, j) = x1;
        end
    end
end

% Save results to a file
save('results', 'triv_0_list', 'triv_1_list', 'N_list', 'g_list');

%% Auxillary functions for fsolve

% Function to compute the maximum real part of the eigenvalues of the Jacobian matrix
% for the trivial_zero_stab stability criterion.
%
% This function calculates the maximum real part of the eigenvalues of the Jacobian
% matrix for the trivial_zero_stab stability criterion based on the given length
% parameter (LI) and number of elements (N).
%
% Parameters:
% LI - Length parameter value
% p  - Parameter vector containing fixed parameters and current g value
% N  - Number of elements
%
% Returns:
% val - Maximum real part of the eigenvalues of the Jacobian matrix
function val = trivial_zero_stab(LI, p, N)
    y = zeros(2*N, 1);       % Initialize state vector
    p(2) = LI;               % Update parameter vector with current LI value
    val = max(real(eig(SL_jac(0, y, p, N)))); % Compute maximum real part of eigenvalues
end

% Function to compute the maximum real part of the eigenvalues of the Jacobian matrix
% for the trivial_one_stab stability criterion.
%
% This function calculates the maximum real part of the eigenvalues of the Jacobian
% matrix for the trivial_one_stab stability criterion based on the given length
% parameter (LI) and number of elements (N).
%
% Parameters:
% LI - Length parameter value
% p  - Parameter vector containing fixed parameters and current g value
% N  - Number of elements
%
% Returns:
% val - Maximum real part of the eigenvalues of the Jacobian matrix
function val = trivial_one_stab(LI, p, N)
    y = [ones(N,1); zeros(N-1, 1)]; % Initialize state vector
    p(2) = LI;                     % Update parameter vector with current LI value
    val = max(real(eig(SL_polar_jac(0, y, p, N)))); % Compute maximum real part of eigenvalues
end
