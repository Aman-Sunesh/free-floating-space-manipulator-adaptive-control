clear;
clc;
close all;

% Project setup
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

cfg = paper_params();


%% =========================================================
% Proof-of-concept simulation settings
%
% We use zero initial angular momentum:
%
%       p0 = 0
%
% Only rotational dynamics are simulated. The visualization fixes
% the spacecraft COM at the origin and omits base translation.
%
% CASE 1:
%   Prescribed joint motion
%   Hbm*phi_dot is generally nonzero
%   -> spacecraft rotates in reaction
%
% CASE 2:
%   Same nominal command projected into the RNS
%   Hbm*phi_dot is approximately zero
%   -> spacecraft remains approximately fixed
%% =========================================================

dt = 0.01;
t_final = 8;

t = (0:dt:t_final).';
N = length(t);

p0 = 0;

theta0 = 0;
phi0 = cfg.qm_0;


% Storage
theta_open = zeros(N,1);
theta_rns  = zeros(N,1);

phi_open = zeros(N,3);
phi_rns  = zeros(N,3);

omega_open = zeros(N,1);
omega_rns  = zeros(N,1);

reaction_open = zeros(N,1);
reaction_rns  = zeros(N,1);

ee_open = zeros(N,2);
ee_rns  = zeros(N,2);

theta_open(1) = theta0;
theta_rns(1)  = theta0;

phi_open(1,:) = phi0.';
phi_rns(1,:)  = phi0.';


% Simulate both cases
for k = 1:N-1
    tk = t(k);

    zeta = [
         0.45 * cos(0.70*tk);
        -0.35 * sin(0.90*tk + 0.40);
         0.40 * cos(0.50*tk + 0.80)
    ];


    %% =====================================================
    % CASE 1
    % Reaction-unaware joint motion
    %% =====================================================
    phi_k = phi_open(k,:).';

    [Hb, Hbm, ~] = momentum_matrices( ...
        phi_k, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);


    % Use nominal command directly
    phi_dot = zeta;


    % Xu et al. (2016), Eq. (2): Angular-momentum conservation:
    %
    % p0 = Hb*omega_b + Hbm*phi_dot
    %
    % therefore
    %
    % omega_b = (p0 - Hbm*phi_dot)/Hb

    omega_b = ...
        (p0 - Hbm*phi_dot) / Hb;


    % Manipulator reaction contribution
    reaction = Hbm * phi_dot;


    % Euler integration
    theta_open(k+1) = ...
        theta_open(k) + dt*omega_b;

    phi_open(k+1,:) = ...
        (phi_k + dt*phi_dot).';


    omega_open(k) = omega_b;
    reaction_open(k) = reaction;


    %% =====================================================
    % CASE 2
    % Reaction-null-space joint motion
    %% =====================================================

    phi_k = phi_rns(k,:).';

    [Hb, Hbm, ~] = momentum_matrices( ...
        phi_k, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);


    % Xu et al. (2016), Eq. (7): Reaction-null-space projector
    T_rns = reaction_null_projector(Hbm);


    % Project the SAME nominal command into the RNS
    phi_dot = T_rns * zeta;


    %% -----------------------------------------------------
    % Keep the projected arm visibly active
    %
    % Scaling does not destroy the RNS property:
    %
    % Hbm*(c*T*zeta)
    %   = c*Hbm*T*zeta
    %   = 0
    %% -----------------------------------------------------

    nominal_speed = norm(zeta);
    projected_speed = norm(phi_dot);

    if projected_speed > 1e-8

        scale = nominal_speed / projected_speed;

        % Avoid excessive amplification
        scale = min(scale, 2.0);

        phi_dot = scale * phi_dot;

    end


    % Same momentum conservation equation
    omega_b = ...
        (p0 - Hbm*phi_dot) / Hb;


    reaction = Hbm * phi_dot;


    % Euler integration
    theta_rns(k+1) = ...
        theta_rns(k) + dt*omega_b;

    phi_rns(k+1,:) = ...
        (phi_k + dt*phi_dot).';


    omega_rns(k) = omega_b;
    reaction_rns(k) = reaction;

end


% Final logging samples
omega_open(end) = omega_open(end-1);
omega_rns(end)  = omega_rns(end-1);

reaction_open(end) = reaction_open(end-1);
reaction_rns(end)  = reaction_rns(end-1);


% Compute end-effector trajectories for trails
for k = 1:N

    joints = joint_positions( ...
        theta_open(k), ...
        phi_open(k,:).', ...
        cfg.L, ...
        cfg.r(1));

    ee_open(k,:) = joints(:,end).';


    joints = joint_positions( ...
        theta_rns(k), ...
        phi_rns(k,:).', ...
        cfg.L, ...
        cfg.r(1));

    ee_rns(k,:) = joints(:,end).';

end


% Visualization
fig = figure( ...
    'Name', 'Reaction Null Space Proof of Concept', ...
    'Color', 'w', ...
    'Position', [100 100 1400 650]);


tl = tiledlayout( ...
    fig, ...
    1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');


% View limits
plot_lim = 5.5;

% Spacecraft drawing dimensions
base_width = 1.6;
base_height = 0.9;


%% =========================================================
% LEFT
% Reaction-unaware motion
%% =========================================================
ax1 = nexttile(tl,1);

hold(ax1,'on');
box(ax1,'on');
axis(ax1,'equal');

xlim(ax1,[-plot_lim plot_lim]);
ylim(ax1,[-plot_lim plot_lim]);

set(ax1, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k');

xlabel(ax1,'x (m)');
ylabel(ax1,'y (m)');

title(ax1, ...
    'Reaction-Unaware Joint Motion', ...
    'FontWeight','bold');


% Spacecraft body
base1 = patch( ...
    ax1, ...
    NaN, NaN, ...
    [0.80 0.80 0.80], ...
    'EdgeColor',[0.25 0.25 0.25], ...
    'LineWidth',1.5);


% Spacecraft COM -> first joint
mount1 = plot( ...
    ax1, ...
    NaN, NaN, ...
    '-', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',3);


% Manipulator
arm1 = plot( ...
    ax1, ...
    NaN, NaN, ...
    '-o', ...
    'Color',[0.85 0.15 0.15], ...
    'MarkerFaceColor','w', ...
    'MarkerSize',7, ...
    'LineWidth',3);


% EE trail
trail1 = plot( ...
    ax1, ...
    NaN, NaN, ...
    ':', ...
    'Color',[0.85 0.15 0.15], ...
    'LineWidth',1.2);


% Spacecraft COM
plot( ...
    ax1, ...
    0, 0, ...
    'o', ...
    'Color',[0.15 0.15 0.15], ...
    'MarkerFaceColor','w', ...
    'MarkerSize',6);


% Live numerical information
info1 = text( ...
    ax1, ...
    0.03, 0.97, '', ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',11, ...
    'BackgroundColor','w');


%% =========================================================
% RIGHT
% Reaction-null-space motion
%% =========================================================

ax2 = nexttile(tl,2);

hold(ax2,'on');
box(ax2,'on');
axis(ax2,'equal');

xlim(ax2,[-plot_lim plot_lim]);
ylim(ax2,[-plot_lim plot_lim]);

set(ax2, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k');

xlabel(ax2,'x (m)');
ylabel(ax2,'y (m)');

title(ax2, ...
    'Reaction-Null-Space Joint Motion', ...
    'FontWeight','bold');


base2 = patch( ...
    ax2, ...
    NaN, NaN, ...
    [0.80 0.80 0.80], ...
    'EdgeColor',[0.25 0.25 0.25], ...
    'LineWidth',1.5);


mount2 = plot( ...
    ax2, ...
    NaN, NaN, ...
    '-', ...
    'Color',[0.35 0.35 0.35], ...
    'LineWidth',3);


arm2 = plot( ...
    ax2, ...
    NaN, NaN, ...
    '-o', ...
    'Color',[0.05 0.35 0.85], ...
    'MarkerFaceColor','w', ...
    'MarkerSize',7, ...
    'LineWidth',3);


trail2 = plot( ...
    ax2, ...
    NaN, NaN, ...
    ':', ...
    'Color',[0.05 0.35 0.85], ...
    'LineWidth',1.2);


plot( ...
    ax2, ...
    0, 0, ...
    'o', ...
    'Color',[0.15 0.15 0.15], ...
    'MarkerFaceColor','w', ...
    'MarkerSize',6);


info2 = text( ...
    ax2, ...
    0.03, 0.97, '', ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',11, ...
    'BackgroundColor','w');


% Animation
frame_stride = 4;

for k = 1:frame_stride:N

    %% -----------------------------------------------------
    % LEFT: reaction-unaware
    %% -----------------------------------------------------

    Rb = [
        cos(theta_open(k)), -sin(theta_open(k));
        sin(theta_open(k)),  cos(theta_open(k))
    ];


    % Spacecraft rectangle in its body frame
    body_local = [
        -base_width/2,  base_width/2, ...
         base_width/2, -base_width/2;

        -base_height/2, -base_height/2, ...
         base_height/2,  base_height/2
    ];


    base_xy = Rb * body_local;


    joints = joint_positions( ...
        theta_open(k), ...
        phi_open(k,:).', ...
        cfg.L, ...
        cfg.r(1));


    set(base1, ...
        'XData',base_xy(1,:), ...
        'YData',base_xy(2,:));


    set(mount1, ...
        'XData',[0 joints(1,1)], ...
        'YData',[0 joints(2,1)]);


    set(arm1, ...
        'XData',joints(1,:), ...
        'YData',joints(2,:));


    set(trail1, ...
        'XData',ee_open(1:k,1), ...
        'YData',ee_open(1:k,2));


    set(info1, ...
        'String',sprintf([ ...
        't = %.2f s\n' ...
        '\\theta_b = %+.3f deg\n' ...
        '\\omega_b = %+.3f deg/s\n' ...
        'H_{bm}\\phi_dot = %+.3f'], ...
        t(k), ...
        rad2deg(theta_open(k)), ...
        rad2deg(omega_open(k)), ...
        reaction_open(k)));


    %% -----------------------------------------------------
    % RIGHT: reaction-null-space
    %% -----------------------------------------------------

    Rb = [
        cos(theta_rns(k)), -sin(theta_rns(k));
        sin(theta_rns(k)),  cos(theta_rns(k))
    ];


    base_xy = Rb * body_local;


    joints = joint_positions( ...
        theta_rns(k), ...
        phi_rns(k,:).', ...
        cfg.L, ...
        cfg.r(1));


    set(base2, ...
        'XData',base_xy(1,:), ...
        'YData',base_xy(2,:));


    set(mount2, ...
        'XData',[0 joints(1,1)], ...
        'YData',[0 joints(2,1)]);


    set(arm2, ...
        'XData',joints(1,:), ...
        'YData',joints(2,:));


    set(trail2, ...
        'XData',ee_rns(1:k,1), ...
        'YData',ee_rns(1:k,2));


    set(info2, ...
        'String',sprintf([ ...
        't = %.2f s\n' ...
        '\\theta_b = %+.6f deg\n' ...
        '\\omega_b = %+.6f deg/s\n' ...
        'H_{bm}\\phi_dot = %+.3e'], ...
        t(k), ...
        rad2deg(theta_rns(k)), ...
        rad2deg(omega_rns(k)), ...
        reaction_rns(k)));


    drawnow;

    pause(frame_stride*dt);

end


% Numerical proof
fprintf('\n========================================\n');
fprintf('REACTION-NULL-SPACE PROOF OF CONCEPT\n');
fprintf('========================================\n');


fprintf('\nMaximum absolute spacecraft attitude:\n');

fprintf('Reaction-unaware:\n');
fprintf('    %.6f deg\n', ...
    max(abs(rad2deg(theta_open))));

fprintf('Reaction-null-space:\n');
fprintf('    %.6e deg\n', ...
    max(abs(rad2deg(theta_rns))));


fprintf('\nMaximum absolute reaction Hbm*phi_dot:\n');

fprintf('Reaction-unaware:\n');
fprintf('    %.6f\n', ...
    max(abs(reaction_open)));

fprintf('Reaction-null-space:\n');
fprintf('    %.6e\n', ...
    max(abs(reaction_rns)));


% Spacecraft attitude comparison
figure( ...
    'Name','RNS Proof - Spacecraft Attitude', ...
    'Color','w');

hold on;
box on;


plot( ...
    t, ...
    rad2deg(theta_open), ...
    '-', ...
    'Color',[0.85 0.15 0.15], ...
    'LineWidth',1.5);


plot( ...
    t, ...
    rad2deg(theta_rns), ...
    '-', ...
    'Color',[0.05 0.35 0.85], ...
    'LineWidth',1.5);


xlabel('Time (s)');
ylabel('Spacecraft attitude (degree)');


lgd = legend( ...
    'reaction-unaware motion', ...
    'reaction-null-space motion', ...
    'Location','best');

set(lgd, ...
    'Color','w', ...
    'TextColor','k');

set(gca, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k');


% Reaction comparison
figure( ...
    'Name','RNS Proof - Reaction Momentum', ...
    'Color','w');

hold on;
box on;


plot( ...
    t, ...
    reaction_open, ...
    '-', ...
    'Color',[0.85 0.15 0.15], ...
    'LineWidth',1.5);


plot( ...
    t, ...
    reaction_rns, ...
    '-', ...
    'Color',[0.05 0.35 0.85], ...
    'LineWidth',1.5);


xlabel('Time (s)');
ylabel('H_{bm}\phi_dot');


lgd = legend( ...
    'reaction-unaware motion', ...
    'reaction-null-space motion', ...
    'Location','best');

set(lgd, ...
    'Color','w', ...
    'TextColor','k');

set(gca, ...
    'Color','w', ...
    'XColor','k', ...
    'YColor','k');


%% =========================================================
% Local helper
%
% Returns:
%
% joints(:,1) = J1
% joints(:,2) = J2
% joints(:,3) = J3
% joints(:,4) = end effector
%% =========================================================

function joints = joint_positions(theta_b, phi, L, r0)

    Rb = [
        cos(theta_b), -sin(theta_b);
        sin(theta_b),  cos(theta_b)
    ];


    % Spacecraft COM -> first revolute joint
    J1 = Rb * [r0; 0];


    % Absolute link orientations
    a1 = theta_b + phi(1);

    a2 = ...
        theta_b ...
        + phi(1) ...
        + phi(2);

    a3 = ...
        theta_b ...
        + phi(1) ...
        + phi(2) ...
        + phi(3);


    % Joint positions
    J2 = ...
        J1 ...
        + L(1) * [cos(a1); sin(a1)];


    J3 = ...
        J2 ...
        + L(2) * [cos(a2); sin(a2)];


    EE = ...
        J3 ...
        + L(3) * [cos(a3); sin(a3)];


    joints = [
        J1, ...
        J2, ...
        J3, ...
        EE
    ];

end