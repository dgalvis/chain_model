% Script to visualize stability regions based on numerical results.
%
% This script loads precomputed stability results from a .mat file and creates
% two types of plots to visualize the stability criteria:
% 
% 1. **Stability Curves vs G Parameter**:
%    - Plots the curves for trivial_zero_stab and trivial_one_stab stability criteria
%      against the parameter `g` for a fixed number of elements.
%    - Highlights regions where the system is considered stable (above or below the curves).
%
% 2. **Stability Points vs Number of OSCILLATORS**:
%    - Plots discrete points for trivial_zero_stab and trivial_one_stab stability criteria
%      against the number of elements `N` for a fixed g parameter.
%    - Highlights regions where the system is considered stable (above or below the curves).
%
% The script also customizes the plots with axis labels, titles, legends, and visual enhancements
% such as shaded regions to indicate stability regions.
%
% Dependencies:
% - The script requires the 'results.mat' file with variables 'triv_0_list', 'triv_1_list',
%   'N_list', and 'g_list'.
%
% Outputs:
% - Two figures are generated:
%   1. A plot showing stability curves as a function of  `g`.
%   2. A scatter plot showing stability points as a function of number of elements `N`.



%% Initialise
clear; clc; close all;
addpath('../functions'); % Add the 'functions' directory to the search path

%% Load the data
% Load the results from the saved .mat file
load('results.mat', 'triv_0_list', 'triv_1_list', 'N_list', 'g_list');

%%
Nidx = 1; % Index for the number of elements
figure(); hold all;

% Plot the stability curves with axes switched
plot(triv_0_list(Nidx, :), g_list, 'k', 'linewidth', 2); % Trivial zero stability curve
plot(triv_1_list(Nidx, :), g_list, 'b', 'linewidth', 2); % Trivial one stability curve

% Identify valid data points
validIdx = ~isnan(triv_0_list(Nidx, :)); 

% Define the upper and lower boundaries for filling regions
one_max = 1;   % Upper limit for the region above triv_1 curve
zero_min = -10; % Lower limit for the region below triv_0 curve

% Fill the region to the right of triv_1 curve (now transposed)
region_above = fill([one_max * ones(size(triv_1_list(Nidx, :))), fliplr(triv_1_list(Nidx, :))], ...
    [g_list, fliplr(g_list)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Fill the region to the left of triv_0 curve (now transposed)
region_below = fill([triv_0_list(Nidx, validIdx), zero_min * ones(size(triv_0_list(Nidx, validIdx)))], ...
    [g_list(validIdx), fliplr(g_list(validIdx))], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Set the x-axis limits (switched from y-axis in original)
xlim([zero_min, one_max]);

% Customize the plot
title(['N = ', num2str(N_list(Nidx))]); % Title with current number of elements
xlabel('\lambda_{(i)}'); % New X-axis label
ylabel('g'); % New Y-axis label
set(gca, 'FontSize', 24, 'LineWidth', 1.5); % Customize axis properties
box on; % Add a box around the plot
%legend([region_above, region_below], {'Stable sync', 'Stable zero'}, 'Location', 'Southeast'); % Add legend

box on;
exportgraphics(gcf, 'N4_trivsync_transposed.eps', 'Resolution', 300);


%%
Nidx = 2; % Index for the number of elements
figure(); hold all;

% Plot the stability curves with axes switched
plot(triv_0_list(Nidx, :), g_list, 'k', 'linewidth', 2); % Trivial zero stability curve
plot(triv_1_list(Nidx, :), g_list, 'b', 'linewidth', 2); % Trivial one stability curve

% Identify valid data points
validIdx = ~isnan(triv_0_list(Nidx, :)); 

% Define the upper and lower boundaries for filling regions
one_max = 1;   % Upper limit for the region above triv_1 curve
zero_min = -10; % Lower limit for the region below triv_0 curve

% Fill the region to the right of triv_1 curve (now transposed)
region_above = fill([one_max * ones(size(triv_1_list(Nidx, :))), fliplr(triv_1_list(Nidx, :))], ...
    [g_list, fliplr(g_list)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Fill the region to the left of triv_0 curve (now transposed)
region_below = fill([triv_0_list(Nidx, validIdx), zero_min * ones(size(triv_0_list(Nidx, validIdx)))], ...
    [g_list(validIdx), fliplr(g_list(validIdx))], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Set the x-axis limits (switched from y-axis in original)
xlim([zero_min, one_max]);

% Customize the plot
title(['N = ', num2str(N_list(Nidx))]); % Title with current number of elements
xlabel('\lambda_{(i)}'); % New X-axis label
ylabel('g'); % New Y-axis label
set(gca, 'FontSize', 24, 'LineWidth', 1.5); % Customize axis properties
box on; % Add a box around the plot
%legend([region_above, region_below], {'Stable sync', 'Stable zero'}, 'Location', 'Southeast'); % Add legend

box on;
exportgraphics(gcf, 'N10_trivsync_transposed.eps', 'Resolution', 300);
