clear;
clc;

% Add project folders
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();

fprintf('========================================\n');
fprintf('TEST: REACTION NULL SPACE\n');
fprintf('========================================\n\n');


%% Current configuration
phi = cfg.qm_0;


%% Get momentum coupling matrix
[~, Hbm, ~] = momentum_matrices( ...
    phi, ...
    cfg.l, ...
    cfg.r, ...
    cfg.m, ...
    cfg.I);


fprintf('Hbm:\n');
disp(Hbm);


% Xu et al. (2016), Eq. (7): Reaction-null-space projector
T = reaction_null_projector(Hbm);

fprintf('RNS projector T:\n');
disp(T);


% Xu et al. (2016), Eq. (7): Hbm * T should be zero
reaction_matrix = Hbm * T;

fprintf('Hbm * T:\n');
disp(reaction_matrix);

err_null = norm(reaction_matrix);

fprintf('Null-space error:\n');
fprintf('    %.3e\n\n', err_null);

assert(err_null < 1e-10, ...
    'T does not project into the reaction null space.');


%% ---------------------------------------------------------
% Test 2:
%
% Projectors should satisfy:
%
% T*T = T
%% ---------------------------------------------------------

err_projector = norm(T*T - T);

fprintf('Projector idempotence error:\n');
fprintf('    %.3e\n\n', err_projector);

assert(err_projector < 1e-10, ...
    'RNS matrix is not a valid projector.');


%% ---------------------------------------------------------
% Test 3:
%
% Xu et al. (2016), Eqs. (6)-(7):
% Arbitrary zeta projected through T must create
% zero angular momentum reaction.
%% ---------------------------------------------------------

zeta = [
     0.40;
    -0.25;
     0.70
    ];

phi_dot_rns = rns_joint_velocity( ...
    Hbm, zeta);

reaction = Hbm * phi_dot_rns;

fprintf('RNS joint velocity:\n');
disp(phi_dot_rns);

fprintf('Hbm * phi_dot_rns:\n');
fprintf('    %.3e\n\n', reaction);

assert(abs(reaction) < 1e-10, ...
    'RNS motion still generates spacecraft reaction.');


fprintf('RNS TEST PASSED.\n');