function [r_C0, r_C1, r_C2, r_C3, r_C] = ...
    com_positions(phi, theta_b, r_C0_in, l, r, m)

    % Joint angles
    q1 = phi(1);
    q2 = phi(2);
    q3 = phi(3);

    % Distance from joint i to the link COM
    l1 = l(1);
    l2 = l(2);
    l3 = l(3);

    % Geometry distances:
    % r0: spacecraft COM -> joint 1
    % r1,r2: link COM -> next joint
    r0 = r(1);
    r1 = r(2);
    r2 = r(3);

    % Absolute link orientations
    alpha1 = theta_b + q1;
    alpha2 = theta_b + q1 + q2;
    alpha3 = theta_b + q1 + q2 + q3;

    % Spacecraft COM
    r_C0 = r_C0_in;

    % Joint 1
    r_J1 = r_C0 + r0 * [cos(theta_b); sin(theta_b)];

    % Link 1 COM
    r_C1 = r_J1 + l1 * [cos(alpha1); sin(alpha1)];

    % Joint 2
    r_J2 = r_C1 + r1 * [cos(alpha1); sin(alpha1)];

    % Link 2 COM
    r_C2 = r_J2 + l2 * [cos(alpha2); sin(alpha2)];

    % Joint 3
    r_J3 = r_C2 + r2 * [cos(alpha2); sin(alpha2)];

    % Link 3 COM
    r_C3 = r_J3 + l3 * [cos(alpha3); sin(alpha3)];

    % Total system COM
    r_C = ( ...
        m(1)*r_C0 + ...
        m(2)*r_C1 + ...
        m(3)*r_C2 + ...
        m(4)*r_C3 ) / sum(m);

end