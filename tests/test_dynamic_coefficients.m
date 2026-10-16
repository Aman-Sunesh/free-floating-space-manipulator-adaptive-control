function test_dynamic_coefficients()
% Numerically reconstruct Xu's ten dynamic coefficients from the
% physical inertia model, without using published values in the fit.
% Validate the reconstructed coefficients on unseen configurations.
% This is numerical identification, not a symbolic derivation.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_root));
    cfg = paper_params();

    % Xu et al. published values -- reference ONLY; never used in fitting.
    ad_xu = [ ...
        11.9952; 12.6126; 4.1234; 4.4982; 6.8544; ...
         2.2409; 69.7260; 35.1779; 15.1638; 3.1686];

    % Generate a set of physically valid, varied configurations and rates.
    % Multiple states are needed because one scalar momentum equation at
    % one instant cannot identify ten unknown coefficients.
    rng(42, 'twister');
    N = 240;
    phi_samples = -2.8 + 5.6*rand(3,N);
    omega_samples = -0.22 + 0.44*rand(1,N);
    phi_dot_samples = -0.4 + 0.8*rand(3,N);

    Y = zeros(N,10);
    p_physical = zeros(N,1);
    for k = 1:N
        phi = phi_samples(:,k);
        omega = omega_samples(k);
        phi_dot = phi_dot_samples(:,k);

        % Xu et al. (2016), Eq. (2): Physical angular momentum
        [Hb, Hbm] = momentum_matrices( ...
            phi, cfg.l, cfg.r, cfg.m, cfg.I);
        p_physical(k) = Hb*omega + Hbm*phi_dot;

        % Xu et al. (2016), Eq. (3): Dynamic regressor
        [~, Yd] = dynamic_regressor(phi, omega, phi_dot);
        Y(k,:) = Yd;
    end

    % Fit the ten coefficients on training states (physical data only).
    Ntrain = 200;
    Ytrain = Y(1:Ntrain,:);
    ptrain = p_physical(1:Ntrain);
    assert(rank(Ytrain) == 10, 'The regressor does not have full column rank.');
    ad_reconstructed = Ytrain \ ptrain;

    % Verify on unseen configurations that the reconstruction generalizes.
    holdout_residual = Y(Ntrain+1:end,:)*ad_reconstructed ...
                     - p_physical(Ntrain+1:end);
    max_holdout = max(abs(holdout_residual));
    coeff_differences = ad_reconstructed - ad_xu;
    max_coeff_diff = max(abs(coeff_differences));

    fprintf('\nDYNAMIC PARAMETER RECONSTRUCTION (10 coefficients)\n');
    fprintf('Training-regressor rank: %d/10\n', rank(Ytrain));
    fprintf('Training-regressor condition number: %.6g\n', cond(Ytrain));
    fprintf(' k   Reconstructed       Xu published       Difference\n');
    for k = 1:10
        fprintf('%2d   %14.9f      %10.4f      %+12.9f\n', ...
            k, ad_reconstructed(k), ad_xu(k), coeff_differences(k));
    end
    fprintf('Maximum coefficient difference: %.9g\n', max_coeff_diff);
    fprintf('Max held-out momentum residual: %.9g kg m^2/s\n', max_holdout);

    % At 4 decimal places, rounding errors are at most ~5e-5.
    assert(max_coeff_diff <= 5.1e-5, ...
        'Reconstructed coefficients differ from Xu beyond 4-place rounding.');
    assert(max_holdout < 1e-9, ...
        'Reconstructed coefficients do not reproduce held-out momentum.');
    fprintf('PASS: all ten coefficients match Xu to published precision.\n');
end
