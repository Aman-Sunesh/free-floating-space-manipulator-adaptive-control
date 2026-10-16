function [dz, aux] = state_derivative(t, z, cfg)
    % Velocity-level simulation of Xu et al.'s adaptive RNS controller.

    % Unpack state
    theta_b = z(1);

    phi = z(2:4);

    x = z(5:6);

    ad_bar_hat = z(7:17);
    ak_bar_hat = z(18:23);


    % Split generalized parameter estimates
    ad_hat = ad_bar_hat(1:10);
    p0_hat = ad_bar_hat(11);

    ak_hat = ak_bar_hat(1:4);
    v0_hat = ak_bar_hat(5:6);


    % Desired end-effector trajectory
    [x_d, x_d_dot] = desired_trajectory(t);


    %% ---------------------------------------------------------
    % Spacecraft attitude error
    %
    % For the planar case, the vector part of the quaternion
    % attitude error becomes sin(theta_error / 2).
    %% ---------------------------------------------------------
    R_b = [
        cos(theta_b), -sin(theta_b);
        sin(theta_b),  cos(theta_b)
    ];
    R_bd = [
        cos(cfg.qbd), -sin(cfg.qbd);
        sin(cfg.qbd),  cos(cfg.qbd)
    ];

    Delta_R_b = R_bd.' * R_b;
    theta_error = atan2(Delta_R_b(2,1), Delta_R_b(1,1));

    % Planar z-component of the error-quaternion vector.
    epsilon_bv = sin(theta_error/2);

    % True plant momentum matrices
    I_actual = cfg.I;

    % Optional terminal-link inertia perturbation (not a full capture model)
    if isfield(cfg, 'capture_time') && ...
       isfield(cfg, 'payload_inertia_scale') && ...
       t >= cfg.capture_time

        I_actual(4) = ...
            cfg.I(4) * cfg.payload_inertia_scale;
    end

    [Hb, Hbm, ~] = momentum_matrices( ...
        phi, cfg.l, cfg.r, cfg.m, I_actual);

    %% ---------------------------------------------------------
    % Xu et al. (2016), Eq. (3): Reconstruct estimated momentum matrices
    %
    % Since:
    %
    %   Yd * ad = Hb*omega_b + Hbm*phi_dot
    %
    % Hb is obtained by setting:
    %
    %   omega_b = 1, phi_dot = 0
    %
    % and each column of Hbm is obtained by setting one
    % joint velocity to 1.
    %% ---------------------------------------------------------

    [~, Yd_Hb] = dynamic_regressor( ...
        phi, 1, zeros(3,1));

    Hb_hat = Yd_Hb * ad_hat;


    Hbm_hat = zeros(1,3);

    for j = 1:3
        e = zeros(3,1);
        e(j) = 1;

        [~, Yd_Hbm] = dynamic_regressor( ...
            phi, 0, e);

        Hbm_hat(j) = Yd_Hbm * ad_hat;
    end


    %% ---------------------------------------------------------
    % Xu et al. (2016), Eq. (5): Reconstruct estimated Jacobians
    %
    % Since:
    %
    %   Yk*ak = Jb*omega_b + Jm*phi_dot
    %
    % we can extract Jb and Jm using unit velocities.
    %% ---------------------------------------------------------

    [~, Yk_Jb] = kinematic_regressor( ...
        theta_b, phi, 1, zeros(3,1));

    Jb_hat = Yk_Jb * ak_hat;


    Jm_hat = zeros(2,3);

    for j = 1:3
        e = zeros(3,1);
        e(j) = 1;

        [~, Yk_Jm] = kinematic_regressor( ...
            theta_b, phi, 0, e);

        Jm_hat(:,j) = Yk_Jm * ak_hat;
    end


    % Xu et al. (2016), Eqs. (17)-(18): Full adaptive joint-velocity controller
    [phi_dot_cmd, zeta] = adaptive_full_controller( ...
        Hb_hat, Hbm_hat, Jm_hat, ...
        p0_hat, v0_hat, ...
        epsilon_bv, ...
        x, x_d, x_d_dot, ...
        cfg.lambda_b, cfg.Lambda_x);


    % Joint-velocity execution: alpha = 1 gives Xu's ideal servo assumption
    if isfield(cfg, 'servo_alpha')
        alpha = cfg.servo_alpha;
    else
        alpha = 1.0;
    end

    phi_dot = alpha * phi_dot_cmd;


    %% ---------------------------------------------------------
    % Actual spacecraft angular velocity
    %
    % Momentum conservation:
    %
    %   Hb*omega_b + Hbm*phi_dot = p0
    %
    % Therefore:
    %
    %   omega_b = (p0 - Hbm*phi_dot) / Hb
    %
    % Xu et al. (2016), Eq. (2): Spacecraft angular velocity from momentum
    %% ---------------------------------------------------------
    omega_b = ...
        (cfg.p0_true - Hbm * phi_dot) / Hb;


    %% ---------------------------------------------------------
    % Actual end-effector velocity
    %
    %   x_dot = Yk*ak + v0
    %
    % Xu et al. (2016), Eq. (5): True end-effector velocity
    %% ---------------------------------------------------------
    [Yk_bar, Yk] = kinematic_regressor( ...
        theta_b, phi, omega_b, phi_dot);

    x_dot = Yk * cfg.ak_true + cfg.v0_true;


    % Xu et al. (2016), Eqs. (12)-(13): Dynamic adaptation error y1
    y1 = ...
        Hb_hat * (omega_b + cfg.lambda_b * epsilon_bv) ...
        + Hbm_hat * (phi_dot - phi_dot_cmd);


    % Xu et al. (2016), Eq. (3): Dynamic regressor at actual joint velocities
    [Yd_bar, ~] = dynamic_regressor( ...
        phi, omega_b, phi_dot);


    % Xu et al. (2016), Eq. (22): Kinematic adaptation error y2
    delta_x = x - x_d;

    delta_x_dot = x_dot - x_d_dot;

    y2 = ...
        Jb_hat * omega_b ...
        + Jm_hat * (phi_dot - phi_dot_cmd) ...
        - (delta_x_dot + cfg.Lambda_x * delta_x);


    % Xu et al. (2016), Eqs. (14) and (23): Parameter adaptation laws
    [ad_bar_hat_dot, ak_bar_hat_dot] = parameter_update( ...
        Yd_bar, y1, ...
        Yk_bar, y2, ...
        cfg.Gamma_d, cfg.Gamma_k);


    % State derivative
    theta_b_dot = omega_b;

    dz = [
        theta_b_dot;
        phi_dot;
        x_dot;
        ad_bar_hat_dot;
        ak_bar_hat_dot
    ];


    % Diagnostics for plotting/debugging
    aux.omega_b       = omega_b;
    aux.phi_dot       = phi_dot;
    aux.phi_dot_cmd   = phi_dot_cmd;

    aux.x_d           = x_d;
    aux.x_d_dot       = x_d_dot;
    aux.delta_x       = delta_x;

    aux.epsilon_bv    = epsilon_bv;

    aux.y1            = y1;
    aux.y2            = y2;

    aux.Hb             = Hb;
    aux.Hbm            = Hbm;

    aux.Hb_hat         = Hb_hat;
    aux.Hbm_hat        = Hbm_hat;
    aux.Jb_hat         = Jb_hat;
    aux.Jm_hat         = Jm_hat;

    aux.p0_hat         = p0_hat;
    aux.v0_hat         = v0_hat;

    aux.zeta           = zeta;

end