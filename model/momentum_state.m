function [p, p_base, p_arm, Hb, Hbm] = ...
    momentum_state(phi, omega_b, phi_dot, l, r, m, I)

% Compute momentum matrices at current configuration
[Hb, Hbm, ~] = momentum_matrices(phi, l, r, m, I);

% Angular momentum contribution from spacecraft angular velocity
p_base = Hb * omega_b;

% Angular momentum contribution from joint velocities
p_arm = Hbm * phi_dot;

% Xu et al. (2016), Eq. (2): Total angular momentum (planar form)
p = p_base + p_arm;

end