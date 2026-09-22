clear;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: SECTION B ATTITUDE CONTROLLER\n');
fprintf('========================================\n\n');


%% Current configuration

theta_b = cfg.qb_0;
phi = cfg.qm_0;

omega_b = cfg.qb_dot_0;
phi_dot = cfg.qm_dot_0;


%% Initial adaptive estimates

ad_bar_hat = cfg.ad_bar_hat_0;

ad_hat = ad_bar_hat(1:10);
p0_hat = ad_bar_hat(11);


%% ---------------------------------------------------------
% Construct estimated Hb
%% ---------------------------------------------------------

[~, Yd_Hb] = dynamic_regressor( ...
    phi, 1, zeros(3,1));

Hb_hat = Yd_Hb * ad_hat;


%% ---------------------------------------------------------
% Construct estimated Hbm
%% ---------------------------------------------------------

Hbm_hat = zeros(1,3);

for j = 1:3

    e = zeros(3,1);
    e(j) = 1;

    [~, Yd_Hbm] = dynamic_regressor( ...
        phi, 0, e);

    Hbm_hat(j) = Yd_Hbm * ad_hat;

end


%% ---------------------------------------------------------
% Attitude error
%% ---------------------------------------------------------

theta_error = theta_b - cfg.qbd;

epsilon_bv = sin(theta_error/2);


%% ---------------------------------------------------------
% Pick arbitrary secondary RNS motion
%
% Section B allows arbitrary zeta.
%% ---------------------------------------------------------

zeta = [
     0.3;
    -0.2;
     0.4
    ];


%% ---------------------------------------------------------
% Section B adaptive attitude controller
%% ---------------------------------------------------------

phi_dot_cmd = adaptive_attitude_controller( ...
    Hb_hat, ...
    Hbm_hat, ...
    p0_hat, ...
    epsilon_bv, ...
    cfg.lambda_b, ...
    zeta);


fprintf('Commanded joint velocity:\n');
disp(phi_dot_cmd);


%% ---------------------------------------------------------
% Test 1: output dimensions
%% ---------------------------------------------------------

assert(all(size(phi_dot_cmd) == [3 1]), ...
    'Section B controller output must be 3x1.');

assert(all(isfinite(phi_dot_cmd)), ...
    'Section B controller produced NaN/Inf.');


%% ---------------------------------------------------------
% Test 2: RNS part creates zero estimated reaction
%% ---------------------------------------------------------

T_hat = reaction_null_projector(Hbm_hat);

phi_rns = T_hat * zeta;

reaction_rns = Hbm_hat * phi_rns;

fprintf('Estimated reaction from RNS component:\n');
fprintf('    %.3e\n\n', reaction_rns);

assert(abs(reaction_rns) < 1e-10, ...
    'Section B RNS component is not reactionless.');


%% ---------------------------------------------------------
% Test 3: attitude component satisfies Eq. (8)
%
% Controller:
%
% phi_dot_cmd =
%     T_hat*zeta
%     + pinv(Hbm_hat) *
%       (p0_hat + Hb_hat*lambda_b*epsilon_bv)
%
% Since Hbm_hat*T_hat*zeta = 0,
%
% Hbm_hat*phi_dot_cmd
%
% should equal
%
% p0_hat + Hb_hat*lambda_b*epsilon_bv
%% ---------------------------------------------------------

lhs = Hbm_hat * phi_dot_cmd;

rhs = ...
    p0_hat ...
    + Hb_hat * cfg.lambda_b * epsilon_bv;

fprintf('Hbm_hat * phi_dot_cmd:\n');
fprintf('    %.10f\n', lhs);

fprintf('Required Eq. (8) value:\n');
fprintf('    %.10f\n\n', rhs);

assert(abs(lhs-rhs) < 1e-10, ...
    'Section B controller does not satisfy Eq. (8).');


%% ---------------------------------------------------------
% Test 4: zeta = 0 case
%
% Pure spacecraft attitude regulation.
%% ---------------------------------------------------------

zeta_zero = zeros(3,1);

phi_dot_attitude_only = adaptive_attitude_controller( ...
    Hb_hat, ...
    Hbm_hat, ...
    p0_hat, ...
    epsilon_bv, ...
    cfg.lambda_b, ...
    zeta_zero);

lhs_zero = Hbm_hat * phi_dot_attitude_only;

assert(abs(lhs_zero-rhs) < 1e-10, ...
    'Attitude-only Section B controller failed.');


fprintf('SECTION B ATTITUDE CONTROLLER TEST PASSED.\n');