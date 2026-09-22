clear;
clc;

% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: KNOWN-PARAMETER CONTROLLER\n');
fprintf('========================================\n\n');


%% Initial state

theta_b = cfg.qb_0;
phi = cfg.qm_0;

omega_b0 = cfg.qb_dot_0;
phi_dot0 = cfg.qm_dot_0;


%% ---------------------------------------------------------
% True momentum
%% ---------------------------------------------------------

[p0, ~, ~, Hb, Hbm] = ...
    momentum_state( ...
        phi, ...
        omega_b0, ...
        phi_dot0, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);


%% ---------------------------------------------------------
% True ak
%% ---------------------------------------------------------

M = sum(cfg.m);

m0 = cfg.m(1);
m1 = cfg.m(2);
m2 = cfg.m(3);
m3 = cfg.m(4);

l1 = cfg.l(1);
l2 = cfg.l(2);
l3 = cfg.l(3);

r0 = cfg.r(1);
r1 = cfg.r(2);
r2 = cfg.r(3);
r3 = cfg.r(4);

L1 = l1 + r1;
L2 = l2 + r2;
L3 = l3 + r3;

ak = [
    m0*r0/M;

    L1 - ...
    (m1*l1 + (m2+m3)*L1)/M;

    L2 - ...
    (m2*l2 + m3*L2)/M;

    L3 - m3*l3/M
    ];


%% ---------------------------------------------------------
% True reduced Jm
%% ---------------------------------------------------------

Jm = zeros(2,3);

for j = 1:3

    e = zeros(3,1);
    e(j) = 1;

    [~, Yk_temp] = kinematic_regressor( ...
        theta_b, phi, 0, e);

    Jm(:,j) = Yk_temp * ak;

end


%% ---------------------------------------------------------
% Initial EE position
%% ---------------------------------------------------------

[~, x_body] = forward_kinematics( ...
    phi, cfg.L);

Rb = [
    cos(theta_b), -sin(theta_b);
    sin(theta_b),  cos(theta_b)
    ];

x = cfg.R_C0_0 ...
    + Rb * ([r0;0] + x_body(1:2));


%% ---------------------------------------------------------
% Compute v0
%% ---------------------------------------------------------

[Jb_geom, Jm_geom] = jacobians( ...
    phi, cfg.L, theta_b, r0);

x_dot0 = ...
    cfg.R_C0_dot_0 ...
    + Jb_geom * omega_b0 ...
    + Jm_geom * phi_dot0;

[~, Yk0] = kinematic_regressor( ...
    theta_b, phi, omega_b0, phi_dot0);

v0 = x_dot0 - Yk0*ak;


%% Desired trajectory

[x_d, x_d_dot] = desired_trajectory(0);


%% Attitude error

theta_error = theta_b - cfg.qbd;

epsilon_bv = sin(theta_error/2);


%% ---------------------------------------------------------
% Controller
%% ---------------------------------------------------------

phi_dot_cmd = known_parameter_controller( ...
    Hb, ...
    Hbm, ...
    Jm, ...
    p0, ...
    v0, ...
    epsilon_bv, ...
    x, ...
    x_d, ...
    x_d_dot, ...
    cfg.lambda_b, ...
    cfg.Lambda_x);


fprintf('Commanded joint velocity:\n');
disp(phi_dot_cmd);


%% Basic validity

assert(all(size(phi_dot_cmd) == [3 1]), ...
    'Controller output should be 3x1.');

assert(all(isfinite(phi_dot_cmd)), ...
    'Controller produced NaN or Inf.');


%% ---------------------------------------------------------
% Check attitude-control momentum condition
%
% Hbm*phi_dot should equal:
%
% p0 + Hb*lambda_b*epsilon_bv
%% ---------------------------------------------------------

lhs = Hbm * phi_dot_cmd;

rhs = ...
    p0 ...
    + Hb * cfg.lambda_b * epsilon_bv;

fprintf('Hbm * phi_dot_cmd:\n');
fprintf('    %.8f\n', lhs);

fprintf('Required attitude momentum term:\n');
fprintf('    %.8f\n\n', rhs);

assert(abs(lhs-rhs) < 1e-8, ...
    'Known controller does not satisfy attitude momentum condition.');


fprintf('KNOWN CONTROLLER TEST PASSED.\n');