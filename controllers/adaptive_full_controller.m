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
    
    % Xu et al. (2016), Eq. (8): Attitude-regulation joint velocity
    phi_att_dot = pinv(Hbm_hat) * (p0_hat + Hb_hat * lambda_b * epsilon_bv);
    
    % Remaining end-effector velocity for RNS tracking
    v_remaining = ...
        - v0_hat ...
        + x_d_dot ...
        - Lambda_x * delta_x ...
        - Jm_hat * phi_att_dot;

    % Xu et al. (2016), Eq. (17): Compute the RNS tracking input
    zeta = pinv(Jm_hat * T_hat) * v_remaining;

    % Xu et al. (2016), Eq. (18): Complete joint-velocity command
    phi_dot_cmd = T_hat * zeta + phi_att_dot;

end