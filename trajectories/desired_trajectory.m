function [x_d, x_d_dot] = desired_trajectory(t)
    % Xu et al. (2016), Section IV: Desired end-effector trajectory
    x_d = [3.7 + 0.3 * cos(pi * t); 
           0.2 + 0.3 * sin(pi * t)];
       
    % Analytical derivative of the desired trajectory
    x_d_dot = [-0.3 * pi * sin(pi * t); 
                0.3 * pi * cos(pi * t)];

end