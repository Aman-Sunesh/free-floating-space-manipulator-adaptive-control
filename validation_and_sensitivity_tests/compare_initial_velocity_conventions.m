% Compare Xu's prescribed initial spacecraft angular velocity with
% the velocity recorded after the first joint-velocity command.
% Only the first plotted sample changes; the simulation is unchanged.

clear;
clc;
close all;

%% Setup
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

%% Run the unchanged, complete adaptive simulation.
cfg = paper_params();
sim = simulate_plant(cfg);
assert(numel(sim.t) > 2, 'Simulation did not return enough samples.');

%% Xu et al. (2016), Eq. (2): Initial angular momentum
[Hb0, Hbm0, ~] = momentum_matrices( ...
    cfg.qm_0, cfg.l, cfg.r, cfg.m, cfg.I);

% PRE-CONTROL physical initial state specified by Xu.
omega_pre = cfg.qb_dot_0;
phi_dot_pre = cfg.qm_dot_0;
p_pre = Hb0 * omega_pre + Hbm0 * phi_dot_pre;

% POST-COMMAND physical initial state recorded by our velocity-level plant.
omega_post = sim.omega_b(1);
phi_dot_post = sim.phi_dot(1,:).';
p_post = Hb0 * omega_post + Hbm0 * phi_dot_post;

% Predict post-command angular velocity INDEPENDENTLY from conserved p0.
omega_post_from_momentum = (p_pre - Hbm0 * phi_dot_post) / Hb0;

% Incorrect hybrid: pre-control spacecraft velocity paired with the
% POST-command joint rates. This hybrid does NOT conserve initial momentum.
p_hybrid = Hb0 * omega_pre + Hbm0 * phi_dot_post;
hybrid_momentum_error = abs(p_hybrid - p_pre);

%% Logging-only alternative: replace ONLY the first point of omega_b.
omega_post_plot = sim.omega_b;
omega_pre_plot = sim.omega_b;
omega_pre_plot(1) = omega_pre;

% The simulated states and every logged velocity after t=0 are unchanged.
max_later_velocity_difference = max(abs( ...
    omega_pre_plot(2:end) - omega_post_plot(2:end)));

% Compare apparent first-sample jumps at the 2-ms next logged sample.
first_jump_post_deg = rad2deg(omega_post_plot(2) - omega_post_plot(1));
first_jump_pre_deg = rad2deg(omega_pre_plot(2) - omega_pre_plot(1));

% If one also (incorrectly) integrated the pre-control angular rate for
% the first time step while using commanded post-control joint rates,
% the first Euler attitude increment would differ by this amount.
first_step_attitude_increment_difference_deg = ...
    rad2deg(cfg.dt * (omega_pre - omega_post));

%% Explicit consistency checks.
assert(abs(p_pre - sim.p0_true) < 1e-8, ...
    'The prescribed initial velocities do not reproduce stored p0.');
assert(abs(p_post - sim.p0_true) < 1e-8, ...
    'The commanded velocities do not preserve stored p0.');
assert(abs(omega_post_from_momentum - omega_post) < 1e-10, ...
    'The post-command omega does not match independent momentum calculation.');
assert(max_later_velocity_difference == 0, ...
    'A logging-only change must not alter velocities after the first sample.');

%% Numerical comparison and report
report = sprintf([ ...
    '\nINITIAL ANGULAR-VELOCITY LOGGING CHECK\n' ...
    'Sample interval: %.3f ms\n' ...
    'Pre-control omega: %.9f deg/s\n' ...
    'Post-command omega: %.9f deg/s\n' ...
    'First-sample difference: %.9f deg/s\n' ...
    'Initial momentum (pre/post): %.12f / %.12f kg m^2/s\n' ...
    'Momentum residuals (pre/post): %.3e / %.3e kg m^2/s\n' ...
    'Invalid hybrid momentum residual: %.3e kg m^2/s\n' ...
    'First jump (post/pre logging): %+.9f / %+.9f deg/s\n' ...
    'Later velocity difference: %.3e deg/s\n' ...
    'Hypothetical first-step attitude difference: %+.9f deg\n'], ...
    1000*cfg.dt, ...
    rad2deg(omega_pre), ...
    rad2deg(omega_post), ...
    rad2deg(omega_post-omega_pre), ...
    p_pre, p_post, ...
    abs(p_pre-sim.p0_true), ...
    abs(p_post-sim.p0_true), ...
    hybrid_momentum_error, ...
    first_jump_post_deg, ...
    first_jump_pre_deg, ...
    rad2deg(max_later_velocity_difference), ...
    first_step_attitude_increment_difference_deg);


fprintf('%s', report);
report_path = fullfile(project_root, 'initial_velocity_comparison.txt');
fid = fopen(report_path, 'w');
assert(fid ~= -1, 'Cannot write the comparison report.');
fprintf(fid, '%s', report);
fclose(fid);

%% Compare the two plots. Same physical trajectory for all t > 0.
fig = figure('Color','w','Name','Xu initial angular velocity conventions');

subplot(2,1,1);
plot(sim.t, rad2deg(omega_post_plot), '-', 'LineWidth', 1.5); hold on;
plot(sim.t, rad2deg(omega_pre_plot), '--', 'LineWidth', 1.2);
plot(sim.t(1), rad2deg(omega_post), 'o', 'MarkerSize', 7);
plot(sim.t(1), rad2deg(omega_pre), 's', 'MarkerSize', 7);
xlim([0, min(0.04,sim.t(end))]); grid on;
ylabel('Angular velocity (deg/s)');
title('First 40 ms: logging convention changes only the first point');
legend('Post-command logging (actual)', 'Prescribed-first logging (display only)', ...
    'Post-command at t=0', 'Pre-control at t=0', 'Location','best');

subplot(2,1,2);
plot(sim.t, rad2deg(omega_post_plot), '-', 'LineWidth', 1.5); hold on;
plot(sim.t, rad2deg(omega_pre_plot), '--', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)'); ylabel('Angular velocity (deg/s)');
title('Full simulation: all subsequent points are identical');
legend('Post-command', 'Prescribed-first (display only)', 'Location','best');

plot_path = fullfile(project_root, 'initial_velocity_comparison.png');
print(fig, plot_path, '-dpng', '-r180');
fprintf('\nSaved report: %s\nSaved figure: %s\n', report_path, plot_path);
