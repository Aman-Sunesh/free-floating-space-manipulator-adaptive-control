function phi_dot_cmd = adaptive_attitude_controller( ...
    Hb_hat, Hbm_hat, p0_hat, epsilon_bv, lambda_b, zeta)

    % Estimated reaction-null-space projector
    T_hat = reaction_null_projector(Hbm_hat);

    % Joint velocity command for spacecraft attitude regulation
    phi_dot_cmd = T_hat * zeta ...
                  + pinv(Hbm_hat) * (p0_hat + Hb_hat * lambda_b * epsilon_bv);

end