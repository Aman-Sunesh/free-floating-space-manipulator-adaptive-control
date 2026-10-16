% Verify numerical convergence of the adaptive controller simulation
% by reducing the integration timestep from 0.004 to 0.00025 s.
% Compare spacecraft angular velocity, attitude, and EE tracking errors
% against the finest-resolution simulation.
clear; clc;
addpath(genpath(pwd));

dts = [0.004, 0.002, 0.001, 0.0005, 0.00025];
sims = cell(size(dts));

for i = 1:length(dts)
    cfg = paper_params();
    cfg.dt = dts(i);

    fprintf('Running dt = %.5f s...\n', cfg.dt);
    sims{i} = simulate_plant(cfg);
end

% Finest result = numerical reference
ref = sims{end};

fprintf('\nDifference from dt = %.5f s reference:\n', dts(end));
fprintf('dt       omega(deg/s)   theta(deg)   max XY error(mm)\n');

for i = 1:length(dts)
    s = sims{i};

    omega_ref = interp1(ref.t, ref.omega_b, s.t);
    theta_ref = interp1(ref.t, ref.theta_b, s.t);

    dx_ref = [
        interp1(ref.t, ref.delta_x(:,1), s.t), ...
        interp1(ref.t, ref.delta_x(:,2), s.t)
    ];

    e_omega = max(abs(rad2deg(s.omega_b - omega_ref)));
    e_theta = max(abs(rad2deg(s.theta_b - theta_ref)));
    e_track = 1000 * max(abs(s.delta_x(:) - dx_ref(:)));

    fprintf('%.5f   %12.6g   %12.6g   %12.6g\n', ...
        dts(i), e_omega, e_theta, e_track);
end