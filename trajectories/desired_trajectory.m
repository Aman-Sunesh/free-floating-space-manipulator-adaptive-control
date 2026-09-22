function [x_d, x_d_dot] = desired_trajectory(t)
    % Desired position trajectory
    x_d = [3.7 + 0.3 * cos(pi * t); 
           0.2 + 0.3 * sin(pi * t)];
       
    % Desired velocity trajectory 
    x_d_dot = [-0.3 * pi * sin(pi * t); 
                0.3 * pi * cos(pi * t)];

end