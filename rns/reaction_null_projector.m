function T = reaction_null_projector(Hbm)
    % Number of manipulator joints
    n = size(Hbm, 2);
    
    % Xu et al. (2016), Eq. (7): Reaction null-space projector
    T = eye(n) - pinv(Hbm) * Hbm;

end