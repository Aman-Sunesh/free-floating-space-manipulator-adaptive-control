function sim = simulate_zero_momentum_plant(cfg)
    % SIMULATE_ZERO_MOMENTUM_PLANT
    %
    % Reproduces the "zero-momenta adaptive controller"
    % used as the comparison controller in Xu et al.
    %
    % IMPORTANT:
    %
    % The TRUE physical plant still has the actual nonzero
    % initial angular momentum p0 and linear momentum / COM
    % velocity v0.
    %
    % The CONTROLLER ignores those quantities:
    %
    %       p0_hat = 0
    %       v0_hat = [0;0]
    %
    % and does not adapt them.
    %
    % Therefore:
    %
    % dynamic estimates  = ad_hat (10 parameters)
    % kinematic estimates = ak_hat (4 parameters)


    %% ---------------------------------------------------------
    % Initial configuration
    %% ---------------------------------------------------------

    theta_b0 = cfg.qb_0;
    phi0 = cfg.qm_0;

    omega_b0 = cfg.qb_dot_0;
    phi_dot0 = cfg.qm_dot_0;


    %% ---------------------------------------------------------
    % Initial EE position from physical geometry
    %% ---------------------------------------------------------

    [~, x_ee_0] = forward_kinematics(phi0, cfg.L);

    Rb0 = [
        cos(theta_b0), -sin(theta_b0);
        sin(theta_b0),  cos(theta_b0)
    ];

    x0 = ...
        cfg.R_C0_0 ...
        + Rb0 * ([cfg.r(1);0] + x_ee_0(1:2));


    %% ---------------------------------------------------------
    % TRUE initial angular momentum
    %
    % Same physical plant as the nonzero-momenta simulation.
    %% ---------------------------------------------------------

    [p0_true, ~, ~] = momentum_state( ...
        phi0, ...
        omega_b0, ...
        phi_dot0, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);

    cfg.p0_true = p0_true;


    %% ---------------------------------------------------------
    % TRUE kinematic parameters
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

    ak1 = m0 * r0 / M;

    ak2 = ...
        L1 ...
        - (m1*l1 + (m2+m3)*L1) / M;

    ak3 = ...
        L2 ...
        - (m2*l2 + m3*L2) / M;

    ak4 = ...
        L3 ...
        - m3*l3 / M;

    cfg.ak_true = [
        ak1;
        ak2;
        ak3;
        ak4
    ];


    %% ---------------------------------------------------------
    % TRUE initial COM velocity term v0
    %% ---------------------------------------------------------

    [Jb_geom0, Jm_geom0] = jacobians( ...
        phi0, ...
        cfg.L, ...
        theta_b0, ...
        r0);

    x_dot0 = ...
        cfg.R_C0_dot_0 ...
        + Jb_geom0 * omega_b0 ...
        + Jm_geom0 * phi_dot0;

    [~, Yk0] = kinematic_regressor( ...
        theta_b0, ...
        phi0, ...
        omega_b0, ...
        phi_dot0);

    cfg.v0_true = ...
        x_dot0 - Yk0 * cfg.ak_true;


    %% ---------------------------------------------------------
    % Initial parameter estimates
    %
    % ZERO-MOMENTA controller:
    %
    % Remove p0_hat from generalized dynamic estimate.
    % Remove v0_hat from generalized kinematic estimate.
    %% ---------------------------------------------------------

    ad_hat0 = cfg.ad_bar_hat_0(1:10);

    ak_hat0 = cfg.ak_bar_hat_0(1:4);


    %% ---------------------------------------------------------
    % State
    %
    % z =
    %
    %   theta_b       1
    %   phi           3
    %   x             2
    %   ad_hat       10
    %   ak_hat        4
    %
    % Total = 20 states
    %% ---------------------------------------------------------

    z0 = [
        theta_b0;
        phi0;
        x0;
        ad_hat0;
        ak_hat0
    ];


    %% ---------------------------------------------------------
    % Time vector
    %% ---------------------------------------------------------

    dt = cfg.dt;

    N = floor(cfg.t_final / dt) + 1;

    t = (0:N-1).' * dt;


    %% ---------------------------------------------------------
    % Storage
    %% ---------------------------------------------------------

    Z = zeros(N,20);
    Z(1,:) = z0.';

    omega_b_log = zeros(N,1);
    phi_dot_log = zeros(N,3);

    x_d_log = zeros(N,2);
    delta_x_log = zeros(N,2);

    y1_log = zeros(N,1);
    y2_log = zeros(N,2);

    Hb_hat_log = zeros(N,1);


    %% ---------------------------------------------------------
    % Simulation loop
    %% ---------------------------------------------------------

    for k = 1:N-1

        zk = Z(k,:).';

        [dz, aux] = zero_momentum_state_derivative( ...
            t(k), ...
            zk, ...
            cfg);

        % Explicit Euler integration
        Z(k+1,:) = (zk + dt*dz).';


        % Diagnostics
        omega_b_log(k) = aux.omega_b;
        phi_dot_log(k,:) = aux.phi_dot.';

        x_d_log(k,:) = aux.x_d.';
        delta_x_log(k,:) = aux.delta_x.';

        y1_log(k) = aux.y1;
        y2_log(k,:) = aux.y2.';

        Hb_hat_log(k) = aux.Hb_hat;

    end


    %% ---------------------------------------------------------
    % Final point
    %% ---------------------------------------------------------

    [~, aux] = zero_momentum_state_derivative( ...
        t(end), ...
        Z(end,:).', ...
        cfg);

    omega_b_log(end) = aux.omega_b;
    phi_dot_log(end,:) = aux.phi_dot.';

    x_d_log(end,:) = aux.x_d.';
    delta_x_log(end,:) = aux.delta_x.';

    y1_log(end) = aux.y1;
    y2_log(end,:) = aux.y2.';

    Hb_hat_log(end) = aux.Hb_hat;


    %% ---------------------------------------------------------
    % Package results
    %% ---------------------------------------------------------

    sim.t = t;
    sim.z = Z;

    sim.theta_b = Z(:,1);

    sim.phi = Z(:,2:4);

    sim.x = Z(:,5:6);

    sim.ad_hat = Z(:,7:16);

    sim.ak_hat = Z(:,17:20);

    sim.omega_b = omega_b_log;

    sim.phi_dot = phi_dot_log;

    sim.x_d = x_d_log;

    sim.delta_x = delta_x_log;

    sim.y1 = y1_log;

    sim.y2 = y2_log;

    sim.Hb_hat = Hb_hat_log;


    % Explicitly show that this controller assumes zero
    % initial momenta.
    sim.p0_hat = zeros(N,1);

    sim.v0_hat = zeros(N,2);


    % Actual physical values
    sim.p0_true = cfg.p0_true;

    sim.v0_true = cfg.v0_true;

    sim.ak_true = cfg.ak_true;

end



%% =============================================================
% ZERO-MOMENTUM CONTROLLER STATE DERIVATIVE
%% =============================================================

function [dz, aux] = zero_momentum_state_derivative(t, z, cfg)


    %% ---------------------------------------------------------
    % Unpack state
    %% ---------------------------------------------------------

    theta_b = z(1);

    phi = z(2:4);

    x = z(5:6);

    ad_hat = z(7:16);

    ak_hat = z(17:20);


    %% ---------------------------------------------------------
    % Desired EE trajectory
    %% ---------------------------------------------------------

    [x_d, x_d_dot] = desired_trajectory(t);


    %% ---------------------------------------------------------
    % Spacecraft attitude error
    %% ---------------------------------------------------------

    theta_error = theta_b - cfg.qbd;

    epsilon_bv = cfg.attitude_error_scale * sin(theta_error/2);


    %% ---------------------------------------------------------
    % TRUE plant momentum matrices
    %% ---------------------------------------------------------

    [Hb, Hbm, ~] = momentum_matrices( ...
        phi, ...
        cfg.l, ...
        cfg.r, ...
        cfg.m, ...
        cfg.I);


    %% ---------------------------------------------------------
    % Estimated Hb
    %% ---------------------------------------------------------

    [~, Yd_Hb] = dynamic_regressor( ...
        phi, ...
        1, ...
        zeros(3,1));

    Hb_hat = Yd_Hb * ad_hat;


    %% ---------------------------------------------------------
    % Estimated Hbm
    %% ---------------------------------------------------------

    Hbm_hat = zeros(1,3);

    for j = 1:3

        e = zeros(3,1);
        e(j) = 1;

        [~, Yd_temp] = dynamic_regressor( ...
            phi, ...
            0, ...
            e);

        Hbm_hat(j) = Yd_temp * ad_hat;

    end


    %% ---------------------------------------------------------
    % Estimated Jb
    %% ---------------------------------------------------------

    [~, Yk_Jb] = kinematic_regressor( ...
        theta_b, ...
        phi, ...
        1, ...
        zeros(3,1));

    Jb_hat = Yk_Jb * ak_hat;


    %% ---------------------------------------------------------
    % Estimated Jm
    %% ---------------------------------------------------------

    Jm_hat = zeros(2,3);

    for j = 1:3

        e = zeros(3,1);
        e(j) = 1;

        [~, Yk_temp] = kinematic_regressor( ...
            theta_b, ...
            phi, ...
            0, ...
            e);

        Jm_hat(:,j) = Yk_temp * ak_hat;

    end


    %% ---------------------------------------------------------
    % ZERO-MOMENTA adaptive controller
    %
    % This is the same structural controller but the effects
    % of p0 and v0 are NOT included.
    %
    % Therefore:
    %
    %       p0_hat = 0
    %       v0_hat = [0;0]
    %% ---------------------------------------------------------

    p0_hat = 0;

    v0_hat = zeros(2,1);

    [phi_dot_cmd, zeta] = adaptive_full_controller( ...
        Hb_hat, ...
        Hbm_hat, ...
        Jm_hat, ...
        p0_hat, ...
        v0_hat, ...
        epsilon_bv, ...
        x, ...
        x_d, ...
        x_d_dot, ...
        cfg.lambda_b, ...
        cfg.Lambda_x);


    %% ---------------------------------------------------------
    % Ideal joint velocity servo
    %% ---------------------------------------------------------

    phi_dot = phi_dot_cmd;


    %% ---------------------------------------------------------
    % TRUE spacecraft angular velocity
    %
    % Actual plant STILL contains p0_true.
    %% ---------------------------------------------------------

    omega_b = ...
        (cfg.p0_true - Hbm*phi_dot) / Hb;


    %% ---------------------------------------------------------
    % TRUE EE velocity
    %
    % Actual plant STILL contains v0_true.
    %% ---------------------------------------------------------

    [~, Yk] = kinematic_regressor( ...
        theta_b, ...
        phi, ...
        omega_b, ...
        phi_dot);

    x_dot = ...
        Yk * cfg.ak_true ...
        + cfg.v0_true;


    %% ---------------------------------------------------------
    % Dynamic adaptation error
    %% ---------------------------------------------------------

    y1 = ...
        Hb_hat * ...
        (omega_b + cfg.lambda_b*epsilon_bv) ...
        + Hbm_hat * ...
        (phi_dot - phi_dot_cmd);


    %% ---------------------------------------------------------
    % Dynamic regressor
    %
    % ZERO-MOMENTA version uses Yd only,
    % NOT [Yd -1].
    %% ---------------------------------------------------------

    [~, Yd] = dynamic_regressor( ...
        phi, ...
        omega_b, ...
        phi_dot);


    %% ---------------------------------------------------------
    % Kinematic adaptation error
    %% ---------------------------------------------------------

    delta_x = x - x_d;

    delta_x_dot = x_dot - x_d_dot;

    y2 = ...
        Jb_hat*omega_b ...
        + Jm_hat*(phi_dot - phi_dot_cmd) ...
        - (delta_x_dot + cfg.Lambda_x*delta_x);


    %% ---------------------------------------------------------
    % ZERO-MOMENTA parameter adaptation
    %
    % Only 10 dynamic parameters and 4 kinematic parameters
    % are estimated.
    %% ---------------------------------------------------------

    Gamma_d_zero = ...
        cfg.Gamma_d(1:10,1:10);

    Gamma_k_zero = ...
        cfg.Gamma_k(1:4,1:4);


    ad_hat_dot = ...
        -Gamma_d_zero ...
        * Yd.' ...
        * y1;


    ak_hat_dot = ...
        -Gamma_k_zero ...
        * Yk.' ...
        * y2;


    %% ---------------------------------------------------------
    % State derivative
    %% ---------------------------------------------------------

    dz = [
        omega_b;
        phi_dot;
        x_dot;
        ad_hat_dot;
        ak_hat_dot
    ];


    %% ---------------------------------------------------------
    % Diagnostics
    %% ---------------------------------------------------------

    aux.omega_b = omega_b;

    aux.phi_dot = phi_dot;

    aux.phi_dot_cmd = phi_dot_cmd;

    aux.x_d = x_d;

    aux.x_d_dot = x_d_dot;

    aux.delta_x = delta_x;

    aux.y1 = y1;

    aux.y2 = y2;

    aux.Hb = Hb;

    aux.Hbm = Hbm;

    aux.Hb_hat = Hb_hat;

    aux.Hbm_hat = Hbm_hat;

    aux.Jb_hat = Jb_hat;

    aux.Jm_hat = Jm_hat;

    aux.zeta = zeta;

end