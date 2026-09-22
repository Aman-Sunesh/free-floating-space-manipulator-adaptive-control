clear;
clc;
close all;

%% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();


%% =========================================================
% 1. Nonzero-momenta controller
%    Planar-scaled attitude error:
%
%    epsilon_bv = 2*sin(theta_error/2)
%% =========================================================

fprintf('Running nonzero-momenta controller with scale = 2...\n');

cfg_scaled = cfg;
cfg_scaled.attitude_error_scale = 2;

sim_nonzero_scaled = simulate_plant(cfg_scaled);


%% =========================================================
% 2. Nonzero-momenta controller
%    Literal quaternion-vector attitude error:
%
%    epsilon_bv = sin(theta_error/2)
%% =========================================================

fprintf('Running nonzero-momenta controller with scale = 1...\n');

cfg_literal = cfg;
cfg_literal.attitude_error_scale = 1;

sim_nonzero_literal = simulate_plant(cfg_literal);


%% =========================================================
% 3. Zero-momenta adaptive controller
%
% Use the scale = 2 convention corresponding to the
% reproduced paper curve.
%% =========================================================

fprintf('Running zero-momenta adaptive controller...\n');

cfg_zero = cfg;
cfg_zero.attitude_error_scale = 2;

sim_zero = simulate_zero_momentum_plant(cfg_zero);


%% =========================================================
% Output directory
%% =========================================================

output_dir = fullfile(project_root, 'results', 'figures');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


%% =========================================================
% Plot settings
%% =========================================================

lw = 1.2;

blue  = [0.0, 0.25, 0.80];
red   = [1.0, 0.0, 0.0];
green = [0.0, 0.55, 0.20];


%% =========================================================
% FIGURE
% Spacecraft attitude convention comparison
%% =========================================================

figure( ...
    'Name', 'Spacecraft Attitude - Error Convention Comparison', ...
    'Color', 'w');

hold on;
box on;

set(gca, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k');


% Zero-momenta adaptive controller
plot( ...
    sim_zero.t, ...
    rad2deg(sim_zero.theta_b), ...
    '--', ...
    'Color', blue, ...
    'LineWidth', lw);


% Nonzero-momenta controller:
% epsilon_bv = 2*sin(theta_error/2)
plot( ...
    sim_nonzero_scaled.t, ...
    rad2deg(sim_nonzero_scaled.theta_b), ...
    '-', ...
    'Color', red, ...
    'LineWidth', lw);


% Nonzero-momenta controller:
% epsilon_bv = sin(theta_error/2)
plot( ...
    sim_nonzero_literal.t, ...
    rad2deg(sim_nonzero_literal.theta_b), ...
    '-', ...
    'Color', green, ...
    'LineWidth', lw);


xlim([0 8]);

% Automatic limits so the literal quaternion case
% is not clipped.
ylim auto;

xticks(0:1:8);
yticks auto;

xlabel('Time (s)');
ylabel('Spacecraft attitude (degree)');


lgd = legend( ...
    'zero-momenta adaptive controller', ...
    'nonzero-momenta: 2sin(\Delta\theta/2)', ...
    'nonzero-momenta: sin(\Delta\theta/2)', ...
    'Location', 'northwest');

set(lgd, ...
    'Color', 'w', ...
    'TextColor', 'k');

set(gca, 'FontSize', 11);


%% =========================================================
% Export
%% =========================================================

exportgraphics( ...
    gcf, ...
    fullfile( ...
        output_dir, ...
        'Fig3_attitude_error_convention_comparison.png'), ...
    'Resolution', 300, ...
    'BackgroundColor', 'white');


%% =========================================================
% Summary
%% =========================================================

fprintf('\n========================================\n');
fprintf('ATTITUDE CONVENTION COMPARISON COMPLETE\n');
fprintf('========================================\n');

fprintf('\nFinal attitude, scale = 2:\n');
fprintf('    %.6f deg\n', ...
    rad2deg(sim_nonzero_scaled.theta_b(end)));

fprintf('\nFinal attitude, scale = 1:\n');
fprintf('    %.6f deg\n', ...
    rad2deg(sim_nonzero_literal.theta_b(end)));

fprintf('\nFigure saved to:\n%s\n', ...
    fullfile( ...
        output_dir, ...
        'Fig3_attitude_error_convention_comparison.png'));