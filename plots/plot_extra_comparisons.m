
clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

out = fullfile(project_root, 'results', 'figures');
if ~exist(out, 'dir'), mkdir(out); end

%% 1. Joint-velocity execution mismatch

% Results from the completed alpha sweep
alpha = [0.80 0.90 0.95 0.98 1.00 1.05 1.10 1.20];
error_mm = [10.21300 4.574442 2.178637 0.858236 ...
    0.136550 1.971709 3.769085 6.932853];

figure('Color','w');
plot(alpha, error_mm, 'o-', 'LineWidth', 1.8, ...
    'MarkerSize', 7);
hold on; grid on; box on;

xline(1, 'k--', 'Ideal execution', ...
    'LabelVerticalAlignment','bottom');

xlabel('Joint-velocity execution scale \alpha');
ylabel('Final EE tracking-error norm (mm)');
title('Sensitivity to Joint-Velocity Execution Mismatch');

exportgraphics(gcf, fullfile(out, ...
    'ServoMismatch_final_error.png'), 'Resolution',300);


%% 2. Figure 3: unscaled vs 2x attitude feedback

cfg = paper_params();

% Original quaternion feedback (scale = 1)
nonzero_1 = simulate_plant(cfg);
zero_1 = simulate_zero_momentum_plant(cfg);

% Double the attitude-feedback gain in the controller and adaptation error
cfg.lambda_b = 2 * cfg.lambda_b;

nonzero_2 = simulate_plant(cfg);
zero_2 = simulate_zero_momentum_plant(cfg);

figure('Color','w');
hold on; grid on; box on;

plot(zero_1.t, rad2deg(zero_1.theta_b), ...
    'b-', 'LineWidth',1.5);

plot(zero_2.t, rad2deg(zero_2.theta_b), ...
    'b--', 'LineWidth',1.5);

plot(nonzero_1.t, rad2deg(nonzero_1.theta_b), ...
    'r-', 'LineWidth',1.5);

plot(nonzero_2.t, rad2deg(nonzero_2.theta_b), ...
    'r--', 'LineWidth',1.5);

xlabel('Time (s)');
ylabel('Spacecraft attitude (deg)');
title('Effect of Attitude-Feedback Scaling');

legend('Zero-momenta: 1x', 'Zero-momenta: 2x', ...
    'Nonzero-momenta: 1x', 'Nonzero-momenta: 2x', ...
    'Location','best');

xlim([0 8]);

exportgraphics(gcf, fullfile(out, ...
    'Fig3_scaling_comparison.png'), 'Resolution',300);

fprintf('\nBoth comparison figures saved to:\n%s\n', out);
