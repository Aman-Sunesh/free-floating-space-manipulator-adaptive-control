clear;
clc;

% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: KINEMATICS\n');
fprintf('========================================\n\n');


%% Test configuration

theta_b = 0.20;

phi = [
     0.60;
    -0.80;
     0.40
    ];

omega_b = -0.07;

phi_dot = [
     0.05;
    -0.03;
     0.08
    ];


%% ---------------------------------------------------------
% 1. Test raw geometric Jacobians using finite differences
%% ---------------------------------------------------------

[Jb, Jm] = jacobians( ...
    phi, cfg.L, theta_b, cfg.r(1));

[~, x_ee_body] = forward_kinematics(phi, cfg.L);

Rb = [
    cos(theta_b), -sin(theta_b);
    sin(theta_b),  cos(theta_b)
    ];

x = Rb * ([cfg.r(1); 0] + x_ee_body(1:2));


% Small timestep
h = 1e-7;

theta_b_next = theta_b + omega_b*h;
phi_next = phi + phi_dot*h;

[~, x_ee_body_next] = forward_kinematics( ...
    phi_next, cfg.L);

Rb_next = [
    cos(theta_b_next), -sin(theta_b_next);
    sin(theta_b_next),  cos(theta_b_next)
    ];

x_next = ...
    Rb_next * ...
    ([cfg.r(1); 0] + x_ee_body_next(1:2));


% Numerical velocity
x_dot_num = (x_next - x) / h;

% Jacobian velocity
x_dot_jac = ...
    Jb * omega_b ...
    + Jm * phi_dot;

jacobian_error = norm(x_dot_num - x_dot_jac);

fprintf('Raw Jacobian finite-difference error:\n');
fprintf('    %.3e\n\n', jacobian_error);

assert(jacobian_error < 1e-5, ...
    'Raw geometric Jacobian test failed.');


%% ---------------------------------------------------------
% 2. Derive true kinematic parameters from physical geometry
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

    L1 ...
    - (m1*l1 + (m2+m3)*L1)/M;

    L2 ...
    - (m2*l2 + m3*L2)/M;

    L3 ...
    - m3*l3/M
    ];

fprintf('Derived ak:\n');
disp(ak.');

fprintf('Xu reported approximately:\n');
disp([0.6277, 1.1550, 1.2600, 1.3542]);


%% ---------------------------------------------------------
% 3. Test Yk using finite differences
%
% Reduced EE position relative to total system COM:
%
% x_rel =
%   ak1*u(alpha0)
% + ak2*u(alpha1)
% + ak3*u(alpha2)
% + ak4*u(alpha3)
%% ---------------------------------------------------------

alpha0 = theta_b;
alpha1 = theta_b + phi(1);
alpha2 = theta_b + phi(1) + phi(2);
alpha3 = theta_b + phi(1) + phi(2) + phi(3);

x_rel = ...
    ak(1) * [cos(alpha0); sin(alpha0)] ...
    + ak(2) * [cos(alpha1); sin(alpha1)] ...
    + ak(3) * [cos(alpha2); sin(alpha2)] ...
    + ak(4) * [cos(alpha3); sin(alpha3)];


theta_next = theta_b + omega_b*h;
phi_next = phi + phi_dot*h;

alpha0_next = theta_next;
alpha1_next = theta_next + phi_next(1);
alpha2_next = theta_next + phi_next(1) + phi_next(2);
alpha3_next = theta_next + phi_next(1) ...
    + phi_next(2) + phi_next(3);

x_rel_next = ...
    ak(1) * [cos(alpha0_next); sin(alpha0_next)] ...
    + ak(2) * [cos(alpha1_next); sin(alpha1_next)] ...
    + ak(3) * [cos(alpha2_next); sin(alpha2_next)] ...
    + ak(4) * [cos(alpha3_next); sin(alpha3_next)];


x_dot_num_reduced = ...
    (x_rel_next - x_rel) / h;

[~, Yk] = kinematic_regressor( ...
    theta_b, phi, omega_b, phi_dot);

x_dot_Yk = Yk * ak;

Yk_error = norm( ...
    x_dot_num_reduced - x_dot_Yk);

fprintf('\nYk finite-difference error:\n');
fprintf('    %.3e\n\n', Yk_error);

assert(Yk_error < 1e-5, ...
    'Kinematic regressor Yk test failed.');


fprintf('KINEMATICS TEST PASSED.\n');