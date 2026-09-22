function [p, p_base, p_arm, Hb, Hbm] = ...
    momentum_state(phi, omega_b, phi_dot, l, r, m, I)

% Compute momentum matrices at current configuration
[Hb, Hbm, ~] = momentum_matrices(phi, l, r, m, I);

% Spacecraft contribution to angular momentum
p_base = Hb * omega_b;

% Manipulator contribution to angular momentum
p_arm = Hbm * phi_dot;

% Total angular momentum
p = p_base + p_arm;

end