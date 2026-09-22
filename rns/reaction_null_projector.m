function T = reaction_null_projector(Hbm)
    % Number of manipulator joints
    n = size(Hbm, 2);
    
    T = eye(n) - pinv(Hbm) * Hbm;

end