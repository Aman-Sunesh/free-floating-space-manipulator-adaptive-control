% Verify that the simulated free-floating spacecraft-manipulator
% conserves total angular momentum throughout the trajectory.
% Compare momentum at every timestep against its initial value.

clear; clc;
addpath(genpath(pwd));

cfg = paper_params();
cfg.dt = 0.002;          % Xu paper setting

sim = simulate_plant(cfg);

% Xu et al. (2016), Eq. (2): Verify angular momentum conservation.
L = zeros(length(sim.t),1);

for k = 1:length(sim.t)

    phi = sim.phi(k,:).';
    omega_b = sim.omega_b(k);
    phi_dot = sim.phi_dot(k,:).';

    L(k) = momentum_state( ...
        phi, ...
        omega_b, ...
        phi_dot, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);
end

fprintf('Initial angular momentum : %.10f\n', L(1));
fprintf('Final angular momentum   : %.10f\n', L(end));
fprintf('Maximum residual         : %.3e\n', max(abs(L - L(1))));