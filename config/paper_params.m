function cfg = paper_params()

% Masses
m0 = 61.2;
m1 = 6.3;
m2 = 5.4;
m3 = 5.1;

% Inertias
I0 = 26.1120;
I1 = 1.0290;
I2 = 0.8820;
I3 = 0.8330;

% Geometry
l0 = 0.80;
l1 = 0.70;
l2 = 0.70;
l3 = 0.70;

r0 = 0.80;
r1 = 0.70;
r2 = 0.70;
r3 = 0.70;

% Desired spacecraft attitude
qbd = 0;

% Planar attitude-error scaling used for figure reproduction
% 2 gives epsilon_bv = 2*sin(theta_error/2) ~ theta_error for small angles
attitude_error_scale = 2;

% Simulation
dt = 0.002;
t_final = 8;

% Initial spacecraft state
R_C0_0 = [0; 0];
R_C0_dot_0 = [0.1; 0.1];

qb_0 = 0;
qb_dot_0 = -0.05;

% Initial joint state
qm_0 = [
    pi/3;
    -2*pi/3;
    pi/3
    ];

qm_dot_0 = [
    0.05;
    -0.01;
    0.09
    ];

% Controller gains
lambda_b = 60;
Lambda_x = 20 * eye(2);

Gamma_d = diag([
    30, 30, ...
    10, 10, 10, 10, 10, 10, 10, 10, ...
    1
    ]);

Gamma_k = diag([
    20, 20, 20, 20, ...
    2, 2
    ]);

% Initial dynamic parameter estimates
ad_bar_hat_0 = [
    30;
    20;
    3;
    3;
    3;
    5;
    100;
    60;
    30;
    2;
    0
    ];

% Initial kinematic parameter estimates
ak_bar_hat_0 = [
    2;
    3;
    3;
    3;
    0;
    0
    ];

% Grouped parameters
m = [m0; m1; m2; m3];
I = [I0; I1; I2; I3];

l = [l1; l2; l3];
r = [r0; r1; r2; r3];

% Full manipulator link lengths
L = [
    l1 + r1;
    l2 + r2;
    l3 + r3
    ];

% Package parameters for simulation
cfg.m = m;
cfg.I = I;

cfg.l = l;
cfg.r = r;
cfg.L = L;

cfg.l0 = l0;

cfg.qbd = qbd;

cfg.dt = dt;
cfg.t_final = t_final;

cfg.R_C0_0 = R_C0_0;
cfg.R_C0_dot_0 = R_C0_dot_0;

cfg.qb_0 = qb_0;
cfg.qb_dot_0 = qb_dot_0;

cfg.qm_0 = qm_0;
cfg.qm_dot_0 = qm_dot_0;

cfg.lambda_b = lambda_b;
cfg.Lambda_x = Lambda_x;

cfg.Gamma_d = Gamma_d;
cfg.Gamma_k = Gamma_k;

cfg.ad_bar_hat_0 = ad_bar_hat_0;
cfg.ak_bar_hat_0 = ak_bar_hat_0;

cfg.attitude_error_scale = attitude_error_scale;

end