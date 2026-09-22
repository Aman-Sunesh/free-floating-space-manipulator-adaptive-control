function [Yk_bar, Yk] = kinematic_regressor( ...
    theta_b, phi, omega_b, phi_dot)

    % Joint angles
    q1 = phi(1);
    q2 = phi(2);
    q3 = phi(3);

    % Joint velocities
    q1_dot = phi_dot(1);
    q2_dot = phi_dot(2);
    q3_dot = phi_dot(3);

    % Absolute orientations
    alpha0 = theta_b;
    alpha1 = theta_b + q1;
    alpha2 = theta_b + q1 + q2;
    alpha3 = theta_b + q1 + q2 + q3;

    % Absolute angular velocities
    w0 = omega_b;
    w1 = omega_b + q1_dot;
    w2 = omega_b + q1_dot + q2_dot;
    w3 = omega_b + q1_dot + q2_dot + q3_dot;

    % Kinematic regressor:
    %
    % x_dot - v0 = Yk * ak
    %
    % ak = [ak1; ak2; ak3; ak4]
    Yk = [
        -w0*sin(alpha0), ...
        -w1*sin(alpha1), ...
        -w2*sin(alpha2), ...
        -w3*sin(alpha3);

         w0*cos(alpha0), ...
         w1*cos(alpha1), ...
         w2*cos(alpha2), ...
         w3*cos(alpha3)
    ];

    % Generalized kinematic regressor:
    %
    % Yk_bar * [ak; v0] = Yk*ak + v0 = x_dot
    Yk_bar = [Yk, eye(2)];

end