clear;
clc;

% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: ADAPTIVE CONTROLLER\n');
fprintf('========================================\n\n');


%% ---------------------------------------------------------
% First check parameter-estimate dimensions
%% ---------------------------------------------------------

assert(length(cfg.ad_bar_hat_0) == 11, ...
    'Dynamic estimate must contain 11 parameters.');

assert(length(cfg.ak_bar_hat_0) == 6, ...
    'Kinematic estimate must contain 6 parameters.');

assert(all(size(cfg.Gamma_d) == [11 11]), ...
    'Gamma_d must be 11x11.');

assert(all(size(cfg.Gamma_k) == [6 6]), ...
    'Gamma_k must be 6x6.');


%% ---------------------------------------------------------
% Test adaptation law dimensions directly
%% ---------------------------------------------------------

phi = cfg.qm_0;
theta_b = cfg.qb_0;

omega_b = cfg.qb_dot_0;
phi_dot = cfg.qm_dot_0;

[Yd_bar, ~] = dynamic_regressor( ...
    phi, omega_b, phi_dot);

[Yk_bar, ~] = kinematic_regressor( ...
    theta_b, phi, omega_b, phi_dot);

% Arbitrary test error signals
y1_test = 0.25;

y2_test = [
     0.10;
    -0.20
    ];

[ad_dot_test, ak_dot_test] = ...
    parameter_update( ...
        Yd_bar, ...
        y1_test, ...
        Yk_bar, ...
        y2_test, ...
        cfg.Gamma_d, ...
        cfg.Gamma_k);


assert(all(size(ad_dot_test) == [11 1]), ...
    'Dynamic adaptation output must be 11x1.');

assert(all(size(ak_dot_test) == [6 1]), ...
    'Kinematic adaptation output must be 6x1.');

assert(all(isfinite(ad_dot_test)), ...
    'Dynamic adaptation produced NaN/Inf.');

assert(all(isfinite(ak_dot_test)), ...
    'Kinematic adaptation produced NaN/Inf.');

fprintf('Adaptation-law dimensions: PASS\n\n');


%% ---------------------------------------------------------
% Run complete simulation
%% ---------------------------------------------------------

fprintf('Running complete adaptive simulation...\n\n');

sim = simulate_plant(cfg);


%% ---------------------------------------------------------
% Basic simulation checks
%% ---------------------------------------------------------

assert(length(sim.t) == 4001, ...
    'Expected 4001 samples for 8 s at dt = 0.002 s.');

assert(all(size(sim.phi,2) == 3), ...
    'Expected three joint trajectories.');

assert(all(size(sim.x,2) == 2), ...
    'Expected a 2D EE trajectory.');

assert(all(isfinite(sim.theta_b)), ...
    'Spacecraft attitude contains NaN/Inf.');

assert(all(isfinite(sim.phi(:))), ...
    'Joint trajectory contains NaN/Inf.');

assert(all(isfinite(sim.x(:))), ...
    'EE trajectory contains NaN/Inf.');

assert(all(isfinite(sim.ad_bar_hat(:))), ...
    'Dynamic estimates contain NaN/Inf.');

assert(all(isfinite(sim.ak_bar_hat(:))), ...
    'Kinematic estimates contain NaN/Inf.');


%% ---------------------------------------------------------
% Validate independently derived plant constants
%% ---------------------------------------------------------

fprintf('Derived p0:\n');
fprintf('    %.8f\n', sim.p0_true);

fprintf('Xu p0:\n');
fprintf('    -1.6467\n\n');

fprintf('Derived v0:\n');
disp(sim.v0_true.');

fprintf('Expected approximately:\n');
disp([0.0988, 0.0943]);

fprintf('Derived ak:\n');
disp(sim.ak_true.');

fprintf('Expected approximately:\n');
disp([0.6277, 1.1550, 1.2600, 1.3542]);


assert(abs(sim.p0_true - (-1.6467)) < 5e-4, ...
    'Derived p0 does not match Xu.');

assert(norm( ...
    sim.v0_true - [0.0988;0.0943]) < 5e-4, ...
    'Derived v0 does not match Xu.');

assert(norm( ...
    sim.ak_true ...
    - [0.6277;1.1550;1.2600;1.3542]) < 5e-4, ...
    'Derived ak does not match Xu.');


%% ---------------------------------------------------------
% Print final tracking diagnostics
%
% Do NOT assert that adaptive parameters equal true values.
% Adaptive control does not require parameter convergence.
%% ---------------------------------------------------------

final_attitude_deg = ...
    rad2deg(sim.theta_b(end));

final_tracking_error = ...
    norm(sim.delta_x(end,:));

fprintf('\nFinal spacecraft attitude:\n');
fprintf('    %.6f deg\n', final_attitude_deg);

fprintf('Final EE tracking-error norm:\n');
fprintf('    %.6e m\n\n', final_tracking_error);


fprintf('ADAPTIVE CONTROLLER TEST PASSED.\n');