function sim = simulate_plant(cfg)
    % SIMULATE_PLANT
    %
    % Runs the velocity-level adaptive FFSM simulation using
    % explicit Euler integration with the paper's sampling period.


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

    x0 = cfg.R_C0_0 ...
       + Rb0 * ([cfg.r(1); 0] + x_ee_0(1:2));


    %% ---------------------------------------------------------
    % Compute TRUE initial angular momentum
    %
    % This should reproduce approximately:
    %
    %       p0 = -1.6467
    %% ---------------------------------------------------------
    [p0_true, ~, ~] = momentum_state( ...
        phi0, ...
        omega_b0, ...
        phi_dot0, ...
        cfg.l, cfg.r, cfg.m, cfg.I);

    cfg.p0_true = p0_true;


    %% ---------------------------------------------------------
    % Derive the TRUE kinematic parameters from physical geometry
    %
    % We derive these rather than inserting Xu's printed values.
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
        - (m1*l1 + (m2 + m3)*L1) / M;

    ak3 = ...
        L2 ...
        - (m2*l2 + m3*L2) / M;

    ak4 = ...
        L3 ...
        - (m3*l3) / M;

    cfg.ak_true = [
        ak1;
        ak2;
        ak3;
        ak4
    ];


    %% ---------------------------------------------------------
    % Compute TRUE v0 from the initial physical motion
    %
    % First compute the actual initial EE velocity from the
    % spacecraft-COM description:
    %
    % x_dot =
    %   R_C0_dot
    % + Jb_geom*omega_b
    % + Jm_geom*phi_dot
    %% ---------------------------------------------------------
    [Jb_geom0, Jm_geom0] = jacobians( ...
        phi0, cfg.L, theta_b0, r0);

    x_dot0 = ...
        cfg.R_C0_dot_0 ...
        + Jb_geom0 * omega_b0 ...
        + Jm_geom0 * phi_dot0;


    % From:
    %
    %   x_dot = Yk*ak + v0
    %
    % therefore:
    %
    %   v0 = x_dot - Yk*ak
    [~, Yk0] = kinematic_regressor( ...
        theta_b0, phi0, omega_b0, phi_dot0);

    cfg.v0_true = ...
        x_dot0 - Yk0 * cfg.ak_true;


    %% ---------------------------------------------------------
    % Print important validation values
    %% ---------------------------------------------------------
    fprintf('Derived initial angular momentum p0:\n');
    fprintf('    %.6f\n\n', cfg.p0_true);

    fprintf('Derived kinematic parameters ak:\n');
    disp(cfg.ak_true.');

    fprintf('Derived initial translational term v0:\n');
    disp(cfg.v0_true.');


    %% ---------------------------------------------------------
    % Initial adaptive parameter estimates
    %% ---------------------------------------------------------
    ad_bar_hat0 = cfg.ad_bar_hat_0;
    ak_bar_hat0 = cfg.ak_bar_hat_0;


    %% ---------------------------------------------------------
    % Full initial state
    %% ---------------------------------------------------------
    z0 = [
        theta_b0;
        phi0;
        x0;
        ad_bar_hat0;
        ak_bar_hat0
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
    nz = length(z0);

    Z = zeros(N, nz);
    Z(1,:) = z0.';

    omega_b_log = zeros(N,1);
    phi_dot_log = zeros(N,3);

    x_d_log = zeros(N,2);
    delta_x_log = zeros(N,2);

    y1_log = zeros(N,1);
    y2_log = zeros(N,2);

    p0_hat_log = zeros(N,1);
    v0_hat_log = zeros(N,2);

    Hb_hat_log = zeros(N,1);


    %% ---------------------------------------------------------
    % Simulation loop
    %% ---------------------------------------------------------
    for k = 1:N-1

        zk = Z(k,:).';

        [dz, aux] = state_derivative(t(k), zk, cfg);

        % Explicit Euler integration
        Z(k+1,:) = (zk + dt * dz).';


        % Store diagnostics
        omega_b_log(k) = aux.omega_b;
        phi_dot_log(k,:) = aux.phi_dot.';

        x_d_log(k,:) = aux.x_d.';
        delta_x_log(k,:) = aux.delta_x.';

        y1_log(k) = aux.y1;
        y2_log(k,:) = aux.y2.';

        p0_hat_log(k) = aux.p0_hat;
        v0_hat_log(k,:) = aux.v0_hat.';

        Hb_hat_log(k) = aux.Hb_hat;

    end


    %% ---------------------------------------------------------
    % Evaluate final point for logging
    %% ---------------------------------------------------------
    [~, aux] = state_derivative(t(end), Z(end,:).', cfg);

    omega_b_log(end) = aux.omega_b;
    phi_dot_log(end,:) = aux.phi_dot.';

    x_d_log(end,:) = aux.x_d.';
    delta_x_log(end,:) = aux.delta_x.';

    y1_log(end) = aux.y1;
    y2_log(end,:) = aux.y2.';

    p0_hat_log(end) = aux.p0_hat;
    v0_hat_log(end,:) = aux.v0_hat.';

    Hb_hat_log(end) = aux.Hb_hat;


    %% ---------------------------------------------------------
    % Package results
    %% ---------------------------------------------------------
    sim.t = t;
    sim.z = Z;

    sim.theta_b = Z(:,1);
    sim.phi = Z(:,2:4);
    sim.x = Z(:,5:6);

    sim.ad_bar_hat = Z(:,7:17);
    sim.ak_bar_hat = Z(:,18:23);

    sim.omega_b = omega_b_log;
    sim.phi_dot = phi_dot_log;

    sim.x_d = x_d_log;
    sim.delta_x = delta_x_log;

    sim.y1 = y1_log;
    sim.y2 = y2_log;

    sim.p0_hat = p0_hat_log;
    sim.v0_hat = v0_hat_log;

    sim.Hb_hat = Hb_hat_log;

    sim.p0_true = cfg.p0_true;
    sim.v0_true = cfg.v0_true;
    sim.ak_true = cfg.ak_true;

end