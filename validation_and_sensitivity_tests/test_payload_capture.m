% Test controller sensitivity to an abrupt terminal-link inertia change.
% This is an idealized perturbation, not a physical payload capture:
% link masses, COM locations, and total momentum remain unchanged.

clear;
clc;
close all;

%% Setup
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

output_dir = fullfile(project_root, 'results', 'inertia_perturbation');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

scales = [1, 2, 5, 10];
perturbation_time = 4.0;
zoom_end = 6.0;

n_cases = numel(scales);
sims = cell(1, n_cases);

colors = lines(n_cases);
labels = arrayfun(@(s) sprintf('I3 scale = %.0fx', s), ...
    scales, 'UniformOutput', false);

fprintf('========================================\n');
fprintf('ABRUPT TERMINAL-INERTIA PERTURBATION\n');
fprintf('Perturbation time: %.2f s\n', perturbation_time);
fprintf('========================================\n\n');

%% Run simulations
for k = 1:n_cases

    cfg = paper_params();

    cfg.capture_time = perturbation_time;
    cfg.payload_inertia_scale = scales(k);

    fprintf('Running terminal inertia scale %.1fx...\n', scales(k));

    sims{k} = simulate_plant(cfg);
end

%% Extract results
t = sims{1}.t(:);
n = numel(t);

attitude_deg = zeros(n, n_cases);
tracking_mm = zeros(n, n_cases);

for k = 1:n_cases

    assert(numel(sims{k}.t) == n && ...
        max(abs(sims{k}.t(:) - t)) < 1e-12, ...
        'Simulation time grids must be identical.');

    attitude_deg(:,k) = ...
        rad2deg(sims{k}.theta_b(:));

    tracking_mm(:,k) = ...
        1000 * vecnorm(sims{k}.delta_x, 2, 2);
end

idx_post = t >= perturbation_time;

%% Figure 1: Full spacecraft attitude
fig1 = figure('Color','w', 'Position',[100 100 950 540]);
hold on; grid on; box on;

h = gobjects(n_cases,1);

for k = 1:n_cases
    h(k) = plot(t, attitude_deg(:,k), ...
        'Color', colors(k,:), 'LineWidth', 1.6);
end

xline(perturbation_time, 'k--', 'Inertia change');

xlabel('Time (s)');
ylabel('Spacecraft attitude (deg)');
title('Spacecraft Attitude Under Abrupt Inertia Perturbation');

lgd = legend(h, labels, 'Location','best');
lgd.AutoUpdate = 'off';
xlim([t(1), t(end)]);

exportgraphics(fig1, ...
    fullfile(output_dir, 'attitude_full.png'), ...
    'Resolution', 300);

%% Figure 2: Attitude difference relative to baseline
% This isolates the effect of the inertia perturbation.

fig2 = figure('Color','w', 'Position',[100 100 950 540]);
hold on; grid on; box on;

h = gobjects(n_cases,1);

for k = 1:n_cases

    delta_attitude = ...
        attitude_deg(:,k) - attitude_deg(:,1);

    h(k) = plot(t, delta_attitude, ...
        'Color', colors(k,:), 'LineWidth', 1.6);
end

xline(perturbation_time, 'k--', 'Inertia change');

xlabel('Time (s)');
ylabel('Attitude difference from baseline (deg)');
title('Spacecraft Attitude Sensitivity to Inertia Change');

lgd = legend(h, labels, 'Location','best');
lgd.AutoUpdate = 'off';
xlim([perturbation_time - 0.1, min(zoom_end,t(end))]);

exportgraphics(fig2, ...
    fullfile(output_dir, 'attitude_difference.png'), ...
    'Resolution', 300);

%% Figure 3: Full end-effector tracking error
fig3 = figure('Color','w', 'Position',[100 100 950 540]);
hold on; grid on; box on;

h = gobjects(n_cases,1);

for k = 1:n_cases
    h(k) = plot(t, tracking_mm(:,k), ...
        'Color', colors(k,:), 'LineWidth', 1.6);
end

xline(perturbation_time, 'k--', 'Inertia change');

xlabel('Time (s)');
ylabel('EE tracking-error norm (mm)');
title('End-Effector Tracking Under Abrupt Inertia Perturbation');

lgd = legend(h, labels, 'Location','best');
lgd.AutoUpdate = 'off';
xlim([t(1), t(end)]);

exportgraphics(fig3, ...
    fullfile(output_dir, 'tracking_full.png'), ...
    'Resolution', 300);

%% Figure 4: Zoomed tracking error
fig4 = figure('Color','w', 'Position',[100 100 950 540]);
hold on; grid on; box on;

h = gobjects(n_cases,1);

for k = 1:n_cases
    h(k) = plot(t, tracking_mm(:,k), ...
        'Color', colors(k,:), 'LineWidth', 1.6);
end

xline(perturbation_time, 'k--', 'Inertia change');

xlabel('Time (s)');
ylabel('EE tracking-error norm (mm)');
title('End-Effector Tracking Near the Inertia Change');

lgd = legend(h, labels, 'Location','best');
lgd.AutoUpdate = 'off';

xlim([perturbation_time - 0.1, min(zoom_end,t(end))]);

exportgraphics(fig4, ...
    fullfile(output_dir, 'tracking_zoom.png'), ...
    'Resolution', 300);

%% Post-perturbation performance metrics

peak_attitude = zeros(n_cases,1);
attitude_difference = zeros(n_cases,1);
rms_tracking = zeros(n_cases,1);
tracking_difference = zeros(n_cases,1);

for k = 1:n_cases

    % Maximum absolute attitude after the perturbation
    peak_attitude(k) = max(abs( ...
        attitude_deg(idx_post,k)));

    % Maximum attitude difference from the unmodified plant
    attitude_difference(k) = max(abs( ...
        attitude_deg(idx_post,k) - ...
        attitude_deg(idx_post,1)));

    % RMS of tracking-error norm, in mm
    rms_tracking(k) = sqrt(mean( ...
        tracking_mm(idx_post,k).^2));

    % Maximum difference between EE error vectors
    % relative to the baseline, in mm
    error_vector_difference = ...
        sims{k}.delta_x(idx_post,:) - ...
        sims{1}.delta_x(idx_post,:);

    tracking_difference(k) = ...
        1000 * max(vecnorm(error_vector_difference,2,2));
end

%% Summary table
summary = table( ...
    scales(:), ...
    peak_attitude, ...
    attitude_difference, ...
    rms_tracking, ...
    tracking_difference, ...
    'VariableNames', { ...
        'InertiaScale', ...
        'PeakAttitude_deg', ...
        'MaxAttitudeDifference_deg', ...
        'RMSTracking_mm', ...
        'MaxTrackingDifference_mm' ...
    });

fprintf('\n========================================\n');
fprintf('POST-PERTURBATION PERFORMANCE SUMMARY\n');
fprintf('========================================\n\n');

disp(summary);

%% Save numerical summary
writetable(summary, ...
    fullfile(output_dir, 'inertia_perturbation_summary.csv'));

fprintf('\nAll figures and summary saved to:\n%s\n', output_dir);
fprintf('\nINERTIA PERTURBATION TEST COMPLETE.\n');
