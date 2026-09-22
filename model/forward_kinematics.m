function [T_0_ee, x_ee_0] = forward_kinematics(phi, L)
    % Link lengths
    a1 = L(1);
    a2 = L(2);
    a3 = L(3);

    % Joint angles
    q1 = phi(1);
    q2 = phi(2);
    q3 = phi(3);

    % End-effector position in frame {0}
    x = (a1 * cos(q1)) + (a2 * cos(q1 + q2)) + (a3 * cos(q1 + q2 + q3));
    y = (a1 * sin(q1)) + (a2 * sin(q1 + q2)) + (a3 * sin(q1 + q2 + q3));
    
    % End-effector orientation
    theta = q1 + q2 + q3;

    % End-effector pose: [x; y; orientation]
    x_ee_0 = [x; y; theta];
    
    % Homogeneous transform from the end-effector frame to frame {0}
    T_0_ee = [
        cos(theta), -sin(theta), x;
        sin(theta),  cos(theta), y;
        0,           0,          1
        ];

end