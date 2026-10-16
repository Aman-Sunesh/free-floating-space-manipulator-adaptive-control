function [Hb, Hbm, Hm] = momentum_matrices(phi, l, r, m, I)
    % Constructs the reduced inertia matrix for the planar
    % free-floating spacecraft + 3R manipulator:
    %
    %       H = [ Hb      Hbm
    %             Hbm'    Hm  ]
    %
    % such that the angular momentum is
    %
    %       p = Hb*omega_b + Hbm*phi_dot

    % Joint angles
    q1 = phi(1);
    q2 = phi(2);
    q3 = phi(3);

    % Link geometry
    l1 = l(1);
    l2 = l(2);
    l3 = l(3);

    r0 = r(1);
    r1 = r(2);
    r2 = r(3);


    % Link direction vectors in the spacecraft frame
    u1 = [
        cos(q1);
        sin(q1)
    ];

    u2 = [
        cos(q1 + q2);
        sin(q1 + q2)
    ];

    u3 = [
        cos(q1 + q2 + q3);
        sin(q1 + q2 + q3)
    ];


    %% ---------------------------------------------------------
    % 90-degree rotation matrix
    %
    % S*u gives derivative of u with respect to its angle:
    %
    %   d/dq [cos(q); sin(q)] = [-sin(q); cos(q)]
    %% ---------------------------------------------------------
    S = [
         0, -1;
         1,  0
    ];

    du1 = S * u1;
    du2 = S * u2;
    du3 = S * u3;

    % Spacecraft COM
    s0 = [0; 0];

    % Link 1 COM
    s1 = [r0; 0] ...
       + l1 * u1;

    % Link 2 COM
    s2 = [r0; 0] ...
       + (l1 + r1) * u1 ...
       + l2 * u2;

    % Link 3 COM
    s3 = [r0; 0] ...
       + (l1 + r1) * u1 ...
       + (l2 + r2) * u2 ...
       + l3 * u3;


    % Body COM position derivatives with respect to joint angles
    Ds0 = zeros(2,3);

    Ds1 = [
        l1 * du1, ...
        zeros(2,1), ...
        zeros(2,1)
    ];

    Ds2 = [
        (l1 + r1)*du1 + l2*du2, ...
        l2*du2, ...
        zeros(2,1)
    ];

    Ds3 = [
        (l1 + r1)*du1 + (l2 + r2)*du2 + l3*du3, ...
        (l2 + r2)*du2 + l3*du3, ...
        l3*du3
    ];


    % Total system COM and its joint-angle derivatives
    M = sum(m);

    sC = ( ...
        m(1)*s0 + ...
        m(2)*s1 + ...
        m(3)*s2 + ...
        m(4)*s3 ) / M;

    % Derivative of total COM with respect to joint angles
    DsC = ( ...
        m(1)*Ds0 + ...
        m(2)*Ds1 + ...
        m(3)*Ds2 + ...
        m(4)*Ds3 ) / M;


    % Body COM positions and derivatives relative to total system COM
    rho0 = s0 - sC;
    rho1 = s1 - sC;
    rho2 = s2 - sC;
    rho3 = s3 - sC;

    Drho0 = Ds0 - DsC;
    Drho1 = Ds1 - DsC;
    Drho2 = Ds2 - DsC;
    Drho3 = Ds3 - DsC;


    %% ---------------------------------------------------------
    % Linear-velocity Jacobians
    %
    % Generalized velocity:
    %
    %   nu = [omega_b;
    %         phi_dot]
    %
    % First column = motion due to spacecraft rotation
    % Last 3 columns = motion due to joint velocities
    %% ---------------------------------------------------------
    Jv0 = [S*rho0, Drho0];
    Jv1 = [S*rho1, Drho1];
    Jv2 = [S*rho2, Drho2];
    Jv3 = [S*rho3, Drho3];


    %% ---------------------------------------------------------
    % Angular-velocity Jacobians
    %
    % omega_0 = omega_b
    % omega_1 = omega_b + q1_dot
    % omega_2 = omega_b + q1_dot + q2_dot
    % omega_3 = omega_b + q1_dot + q2_dot + q3_dot
    %% ---------------------------------------------------------
    Jw0 = [1, 0, 0, 0];
    Jw1 = [1, 1, 0, 0];
    Jw2 = [1, 1, 1, 0];
    Jw3 = [1, 1, 1, 1];


    %% ---------------------------------------------------------
    % Full reduced inertia matrix
    %
    % Each body contributes:
    %
    %   m_i * Jv_i' * Jv_i
    % + I_i * Jw_i' * Jw_i
    %% ---------------------------------------------------------
    H = ...
        m(1) * (Jv0.' * Jv0) + I(1) * (Jw0.' * Jw0) + ...
        m(2) * (Jv1.' * Jv1) + I(2) * (Jw1.' * Jw1) + ...
        m(3) * (Jv2.' * Jv2) + I(3) * (Jw2.' * Jw2) + ...
        m(4) * (Jv3.' * Jv3) + I(4) * (Jw3.' * Jw3);

    % Numerical cleanup: inertia matrix should be symmetric
    H = 0.5 * (H + H.');


    % Xu et al. (2016), Eqs. (1)-(2): Inertia blocks for angular momentum
    Hb  = H(1,1);
    Hbm = H(1,2:4);
    Hm  = H(2:4,2:4);

end