function phi_dot_cmd = known_parameter_controller( ...
    Hb, Hbm, Jm, p0, v0, ...
    attitude_error, ...
    x, x_d, x_d_dot, ...
    lambda_b, Lambda_x)

    % End-effector position error
    delta_x = x - x_d;

    % Reaction null-space projector
    T = reaction_null_projector(Hbm)

    % Joint velocity used for spacecraft attitude regulation
    phi_att_dot = pinv(Hbm) * (p0 + Hb * lambda_b * attitude_error)

    % Remaining EE velocity that must be generated through the RNS
    v_remaining = ...
        -v0 ...
        + x_d_dot ...
        - Lambda_x * delta_x ...
        - Jm * phi_att_dot;

    % Choose zeta so the RNS motion generates the remaining EE velocity
    zeta = pinv(Jm * T) * v_remaining;

    % Total commanded joint velocity
    phi_dot_cmd = phi_att_dot + T * zeta;
    
end
