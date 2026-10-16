function phi_dot_cmd = known_parameter_controller( ...
    Hb, Hbm, Jm, p0, v0, ...
    attitude_error, ...
    x, x_d, x_d_dot, ...
    lambda_b, Lambda_x)

    % End-effector position error
    delta_x = x - x_d;

    % Reaction null-space projector
    T = reaction_null_projector(Hbm);

    % Xu et al. (2016), Eq. (8): Attitude-regulation joint velocity
    phi_att_dot = pinv(Hbm) * (p0 + Hb * lambda_b * attitude_error);

    % Remaining end-effector velocity for RNS tracking
    v_remaining = ...
        -v0 ...
        + x_d_dot ...
        - Lambda_x * delta_x ...
        - Jm * phi_att_dot;

    % Xu et al. (2016), Eq. (17): Compute the RNS tracking input
    zeta = pinv(Jm * T) * v_remaining;

    % Xu et al. (2016), Eq. (18): Complete joint-velocity command
    phi_dot_cmd = phi_att_dot + T * zeta;
    
end
