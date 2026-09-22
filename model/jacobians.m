function [Jb, Jm] = jacobians(phi, L, theta_b, r0)
    % End-effector pose relative to joint 1
    [~, x_ee_0] = forward_kinematics(phi, L);

    % Spacecraft rotation from body frame to inertial frame
    R_b = [
        cos(theta_b), -sin(theta_b);
        sin(theta_b),  cos(theta_b)
        ];

    % Derivative of spacecraft rotation with respect to theta_b
    dRb_thetab = [
        -sin(theta_b), -cos(theta_b);
        cos(theta_b), -sin(theta_b)
        ];

    % Vector from spacecraft COM to EE, expressed in body frame
    r_body  = [r0; 0] + x_ee_0(1:2);

    % Base rotation contribution
    Jb = dRb_thetab * r_body; 

    % Manipulator Jacobian contribution in the inertial frame
    Jm = R_b * jacobian_arm(phi, L);

end