% This file runs the main function and does experiments with different
% problem sizes and # of directions.
%
exp_sizes = [12 14 16 18]; % Test the problems with size 2^12, 2^14, 2^16 and 2^18
num_directions = [4 6 8]; % Test different # of directions
npol = 20; % # of poles
n = 2^18; % Calculate the relative error of the original problem with
                % different # of directions on the problem with size 2^18

relerr_ori_tab = zeros(length(exp_sizes), length(num_directions));
relerr_sub_tab = zeros(length(num_directions), npol);

% Loop over different number of directions and problem sizes
for col = 1:length(num_directions)
    for row = 1:length(exp_sizes)
        n = 2^exp_sizes(row);
        [relerr_orig, relerr_sub, ~] = main(n, [], num_directions(col));
        relerr_ori_tab(row, col) = relerr_orig;
    end
end

% Loop over different number of directions for large problem size
for row = 1:length(num_directions)
    [relerr_orig, relerr_sub] = main(n, [], num_directions(row));
    relerr_sub_tab(row, :) = relerr_sub;
end

% Plotting
hold on
plot_symbols = {'-o', '-s', '-*'};
for i = 1:size(relerr_sub_tab, 1)
    plot(1:npol, log10(relerr_sub_tab(i, :)), plot_symbols{i});
end
xlabel('$i$', 'Interpreter', 'latex')
ylabel('$\log_{10}((u_h(t_i) - u_h(t_i)^{rbm})/u_h(t_i)) $', 'Interpreter', 'latex')
legend('m = 4', 'm = 6', 'm = 8')
hold off
