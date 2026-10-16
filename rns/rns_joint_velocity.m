function phi_dot_rns = rns_joint_velocity(Hbm, zeta)
    % Compute the reaction null-space projector
    T = reaction_null_projector(Hbm);

    % Xu et al. (2016), Eq. (7): Project joint velocity into the RNS
    phi_dot_rns = T * zeta;

end
