function [dz, aux] = state_derivative(t, z, cfg)
    % STATE_DERIVATIVE
    %
    % Velocity-level simulation of Xu et al.'s adaptive controller.
    %
    % State:
    %
    % z = [
    %       theta_b;          % 1
    %       phi;              % 3
    %       x;                % 2
    %       ad_bar_hat;       % 11 = [a_d_hat; p0_hat]
    %       ak_bar_hat        % 6  = [a_k_hat; v0_hat]
    %     ]


    %% ---------------------------------------------------------
    % Unpack state
    %% ---------------------------------------------------------
    theta_b = z(1);

    phi = z(2:4);

    x = z(5:6);

    ad_bar_hat = z(7:17);
    ak_bar_hat = z(18:23);


    %% ---------------------------------------------------------
    % Split generalized parameter estimates
    %% ---------------------------------------------------------
    ad_hat = ad_bar_hat(1:10);
    p0_hat = ad_bar_hat(11);

    ak_hat = ak_bar_hat(1:4);
    v0_hat = ak_bar_hat(5:6);


    %% ---------------------------------------------------------
    % Desired EE trajectory
    %% ---------------------------------------------------------
    [x_d, x_d_dot] = desired_trajectory(t);


    %% ---------------------------------------------------------
    % Spacecraft attitude error
    %
    % For the planar case, the vector part of the quaternion
    % attitude error becomes sin(theta_error / 2).
    %% ---------------------------------------------------------
    theta_error = theta_b - cfg.qbd;

    epsilon_bv = cfg.attitude_error_scale * sin(theta_error/2);


    %% ---------------------------------------------------------
    % TRUE plant momentum matrices
    %% ---------------------------------------------------------
    [Hb, Hbm, ~] = momentum_matrices( ...
        phi, cfg.l, cfg.r, cfg.m, cfg.I);


    %% ---------------------------------------------------------
    % ESTIMATED momentum matrices
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
    % ESTIMATED kinematic Jacobians
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


    %% ---------------------------------------------------------
    % Xu Eq. (17)-(18):
    % full adaptive joint-velocity controller
    %% ---------------------------------------------------------
    [phi_dot_cmd, zeta] = adaptive_full_controller( ...
        Hb_hat, Hbm_hat, Jm_hat, ...
        p0_hat, v0_hat, ...
        epsilon_bv, ...
        x, x_d, x_d_dot, ...
        cfg.lambda_b, cfg.Lambda_x);


    %% ---------------------------------------------------------
    % Ideal joint-velocity servo assumption from the paper
    %
    %       phi_dot == phi_dot_cmd
    %% ---------------------------------------------------------
    phi_dot = phi_dot_cmd;


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
    %% ---------------------------------------------------------
    omega_b = ...
        (cfg.p0_true - Hbm * phi_dot) / Hb;


    %% ---------------------------------------------------------
    % Actual end-effector velocity
    %
    % Xu Eq. (5):
    %
    %   x_dot = Yk*ak + v0
    %% ---------------------------------------------------------
    [Yk_bar, Yk] = kinematic_regressor( ...
        theta_b, phi, omega_b, phi_dot);

    x_dot = Yk * cfg.ak_true + cfg.v0_true;


    %% ---------------------------------------------------------
    % Dynamic adaptation error y1
    %
    % Xu Eq. (12)
    %% ---------------------------------------------------------
    y1 = ...
        Hb_hat * (omega_b + cfg.lambda_b * epsilon_bv) ...
        + Hbm_hat * (phi_dot - phi_dot_cmd);


    %% ---------------------------------------------------------
    % Dynamic regressor at the ACTUAL current motion
    %% ---------------------------------------------------------
    [Yd_bar, ~] = dynamic_regressor( ...
        phi, omega_b, phi_dot);


    %% ---------------------------------------------------------
    % Kinematic adaptation error y2
    %
    % Xu Eq. (22):
    %
    % y2 =
    %   Jb_hat*omega_b
    % + Jm_hat*(phi_dot - phi_dot_cmd)
    % - (delta_x_dot + Lambda_x*delta_x)
    %% ---------------------------------------------------------
    delta_x = x - x_d;

    delta_x_dot = x_dot - x_d_dot;

    y2 = ...
        Jb_hat * omega_b ...
        + Jm_hat * (phi_dot - phi_dot_cmd) ...
        - (delta_x_dot + cfg.Lambda_x * delta_x);


    %% ---------------------------------------------------------
    % Parameter adaptation laws
    %
    % Xu Eq. (14) and Eq. (23)
    %% ---------------------------------------------------------
    [ad_bar_hat_dot, ak_bar_hat_dot] = parameter_update( ...
        Yd_bar, y1, ...
        Yk_bar, y2, ...
        cfg.Gamma_d, cfg.Gamma_k);


    %% ---------------------------------------------------------
    % State derivative
    %% ---------------------------------------------------------
    theta_b_dot = omega_b;

    dz = [
        theta_b_dot;
        phi_dot;
        x_dot;
        ad_bar_hat_dot;
        ak_bar_hat_dot
    ];


    %% ---------------------------------------------------------
    % Diagnostics for plotting/debugging
    %% ---------------------------------------------------------
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