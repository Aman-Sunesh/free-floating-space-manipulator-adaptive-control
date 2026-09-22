clear;
clc;

% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: MOMENTUM\n');
fprintf('========================================\n\n');


%% Initial conditions from Xu

phi = cfg.qm_0;

omega_b = cfg.qb_dot_0;

phi_dot = cfg.qm_dot_0;


%% ---------------------------------------------------------
% 1. Compute physical momentum
%% ---------------------------------------------------------

[p, p_base, p_arm, Hb, Hbm] = ...
    momentum_state( ...
        phi, ...
        omega_b, ...
        phi_dot, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);


fprintf('Hb:\n');
disp(Hb);

fprintf('Hbm:\n');
disp(Hbm);

fprintf('Base contribution:\n');
disp(p_base);

fprintf('Manipulator contribution:\n');
disp(p_arm);

fprintf('Total angular momentum:\n');
fprintf('    %.8f\n\n', p);


%% ---------------------------------------------------------
% Compare against Xu
%% ---------------------------------------------------------

p0_xu = -1.6467;

fprintf('Xu reported p0:\n');
fprintf('    %.4f\n\n', p0_xu);

fprintf('Absolute error:\n');
fprintf('    %.3e\n\n', abs(p - p0_xu));

assert(abs(p - p0_xu) < 5e-4, ...
    'Physical momentum does not reproduce Xu p0.');


%% ---------------------------------------------------------
% 2. Test dynamic regressor
%% ---------------------------------------------------------

ad_xu = [
    11.9952;
    12.6126;
    4.1234;
    4.4982;
    6.8544;
    2.2409;
    69.7260;
    35.1779;
    15.1638;
    3.1686
    ];

[~, Yd] = dynamic_regressor( ...
    phi, omega_b, phi_dot);

p_regressor = Yd * ad_xu;

fprintf('Momentum from Yd * ad:\n');
fprintf('    %.8f\n\n', p_regressor);

fprintf('Physical-model vs regressor difference:\n');
fprintf('    %.3e\n\n', abs(p - p_regressor));

assert(abs(p - p_regressor) < 1e-3, ...
    'Dynamic regressor does not match momentum model.');


%% ---------------------------------------------------------
% 3. Basic inertia sanity checks
%% ---------------------------------------------------------

assert(Hb > 0, ...
    'Hb should be positive.');

assert(all(isfinite(Hbm)), ...
    'Hbm contains invalid values.');


fprintf('MOMENTUM TEST PASSED.\n');