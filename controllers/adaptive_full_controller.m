function [phi_dot_cmd, zeta] = adaptive_full_controller( ...
    Hb_hat, Hbm_hat, Jm_hat, ...
    p0_hat, v0_hat, ...
    epsilon_bv, ...
    x, x_d, x_d_dot, ...
    lambda_b, Lambda_x)

    % End-effector position error
    delta_x = x - x_d;
    
    % Estimated reaction null-space projector
    T_hat = reaction_null_projector(Hbm_hat);
    
    % Joint velocity used for spacecraft attitude regulation
    phi_att_dot = pinv(Hbm_hat) * (p0_hat + Hb_hat * lambda_b * epsilon_bv);
    
    % Remaining EE velocity that must be generated through the RNS
    v_remaining = ...
        - v0_hat ...
        + x_d_dot ...
        - Lambda_x * delta_x ...
        - Jm_hat * phi_att_dot;

    % Choose zeta so the RNS motion generates the remaining EE velocity
    zeta = pinv(Jm_hat * T_hat) * v_remaining;

    % Total commanded joint velocity
    phi_dot_cmd = T_hat * zeta + phi_att_dot;

end