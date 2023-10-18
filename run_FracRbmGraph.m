% This script runs the FracRbmGraph function and conducts experiments for 
% different problem sizes and number of directions.

clear; clc;close;

% Global variables for reuse across function calls.
global uk_exact_all;
global P_all; 

uk_exact_all = [];
P_all = [];

% Define the problem sizes and the number of directions to be experimented with.
exp_sizes = [12 14 16 18];  % These correspond to sizes: 2^12, 2^14, 2^16, and 2^18.
num_directions = [4 6 8];
npol = 20;  % Number of poles.

% Initialize tables to store relative errors.
relerr_ori_tab = zeros(length(exp_sizes), length(num_directions)); % Relative error of the original problem
relerr_sub_tab = zeros(length(num_directions), npol); % Relative error of each subproblem

% Experiment with different numbers of directions and problem sizes.
% It's more efficient to start with the largest number of directions first 
% due to the way P is reused.
num_directions_ = sort(num_directions, 'descend'); 
for i = 1:length(exp_sizes)
    n = 2^exp_sizes(i);
    uk_exact_all = [];
    P_all = [];
    for j = 1:length(num_directions)
        relerr_ori = FracRbmGraph(n, [], num_directions_(j));
        relerr_ori_tab(i, end+1-j) = relerr_ori;
    end
end

% Experiment with different numbers of directions using 2 preconditioners 
% and the largest problem size (n = 2^18).
for i = 1:length(num_directions)
    [relerr_orig, relerr_sub] = FracRbmGraph(n, [], num_directions(i), 2);
    relerr_sub_tab(i, :) = relerr_sub;
end

% Plot the relative errors for each subproblem when n = 2^18.
hold on;
plot_symbols = {'-o', '-s', '-*'};
for i = 1:size(relerr_sub_tab, 1)
    plot(1:npol, log10(relerr_sub_tab(i, :)), plot_symbols{i}, 'LineWidth', 1);
end
