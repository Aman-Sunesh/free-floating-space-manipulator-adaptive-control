function [ad_bar_hat_dot, ak_bar_hat_dot] = parameter_update( ...
    Yd_bar, y1, ...
    Yk_bar, y2, ...
    Gamma_d, Gamma_k)

    % Dynamic parameter adaptation law
    ad_bar_hat_dot = -Gamma_d * Yd_bar.' * y1;
    
    % Kinematic parameter adaptation law
    ak_bar_hat_dot = -Gamma_k * Yk_bar.' * y2;

end