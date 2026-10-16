% Check angular-velocity and attitude consistency in Xu's reproduction.
% Verify numerical integration of spacecraft attitude and compare
% published angular velocity against published attitude.
% Digitized screenshots provide approximate numerical references.

addpath(genpath(pwd));
close all;
cfg = paper_params();

fprintf('\n========== A. OUR SIMULATION: KINEMATIC IDENTITY ==========\n');
S = {simulate_plant(cfg), simulate_zero_momentum_plant(cfg)};
names = {'Nonzero-momenta','Zero-momenta'};
for k = 1:2
    s = S{k};
    h = diff(s.t);
    theta_rate = diff(s.theta_b) ./ h;       % actual time step difference
    omega = s.omega_b(1:end-1);              % Euler uses the left endpoint
    mask = s.t(1:end-1) >= 0.20;              % avoid first-command jump for RMS
    error_full = rad2deg(theta_rate - omega);
    error_half = rad2deg(theta_rate / 2 - omega);
    theta_from_omega = s.theta_b(1) + [0; cumsum(h .* omega)];
    fprintf('\n%s:\n',names{k});
    fprintf('  Max |dtheta/dt - omega| (Euler): %.3e deg/s\n',max(abs(error_full)));
    fprintf('  Max |theta - sum(omega*dt)|:      %.3e deg\n',max(abs(rad2deg(s.theta_b - theta_from_omega))));
    fprintf('  RMS rate residual, true angle:   %.3e deg/s\n',sqrt(mean(error_full(mask).^2)));
    fprintf('  RMS rate residual, HALF angle:   %.6g deg/s\n',sqrt(mean(error_half(mask).^2)));
    fprintf('  Minimum attitude (original):    %.6f deg\n',min(rad2deg(s.theta_b)));
    fprintf('  Minimum attitude (halved):      %.6f deg\n',min(rad2deg(s.theta_b))/2);
    fprintf('  First logged angular velocity:  %.6f deg/s\n',rad2deg(s.omega_b(1)));
end
fprintf('\nPrescribed pre-control omega: %.6f deg/s\n',rad2deg(cfg.qb_dot_0));
fprintf('Note: dtheta/dt = omega is enforced in our state_derivative.\n');

% If original figure images are not available, stop after verifying the code.
f2 = 'actual_fig2.png'; f3 = 'actual_fig3.png';
if ~isfile(f2) || ~isfile(f3)
    fprintf('\n========== B. PUBLISHED FIGURES: NOT YET TESTED ==========\n');
    fprintf('Place %s and %s in this directory, then rerun:\n%s\n',f2,f3,pwd);
    fprintf('The original figure images are needed to assess Xu''s own consistency.\n');
else
    fprintf('\n========== B. DIGITIZE THE PUBLISHED RED CURVES ==========\n');
    % Calibrations for the exact uploaded actual_fig2.png and actual_fig3.png.
    % Pixel axes: [xLeft yTop xRight yBottom]; values: [yTop yBottom].
    [tw,paper_omega,coverage2] = xuDigitizeRed(f2,[136 32 1061 763],[1 -3.5], [141 250],[1118 855]);
    [ta,paper_theta,coverage3] = xuDigitizeRed(f3,[199 23 1361 939],[0.015 -0.030],[170 750],[1403 1032]);
    fprintf('Usable red trace columns: Fig2 %.1f%%, Fig3 %.1f%%\n',100*coverage2,100*coverage3);
    if coverage2 < 0.80 || coverage3 < 0.80
        error('Red trace extraction poor: check that you used the exact actual_fig2/3 PNG images.');
    end

    % Compare only t=0.25..7.75: avoids Fig2 initial vertical spike,
    % annotated red/orange circle, and plot boundary effects.
    t_eval = (0.25:0.005:7.75).';
    w_pub = movmean(interp1(tw,paper_omega,t_eval,'linear'),11);
    th_pub = movmean(interp1(ta,paper_theta,t_eval,'linear'),11);
    integral_w = cumtrapz(t_eval,w_pub);

    % Test the published figure pair: theta(t) = c*int(omega) + offset.
    % Allow a small linear drift nuisance term to compensate for baseline
    % bias of pixel readings from a figure (never modify the source images).
    X = [integral_w, ones(size(t_eval)), t_eval-t_eval(1)];
    fit = X \ th_pub;
    fitted_scale = fit(1);
    drift = fit(3);
    N = [ones(size(t_eval)), t_eval-t_eval(1)];
    d1 = N \ (th_pub - integral_w);
    d05 = N \ (th_pub - 0.5*integral_w);
    rmse1 = sqrt(mean((th_pub - (integral_w + N*d1)).^2));
    rmse05 = sqrt(mean((th_pub - (0.5*integral_w + N*d05)).^2));
    rmsebest = sqrt(mean((th_pub - X*fit).^2));
    fprintf('Best ratio published theta change / integral(omega): %.3f\n',fitted_scale);
    fprintf('Estimated trace baseline drift:                  %.5g deg/s\n',drift);
    fprintf('RMSE if consistent scale=1:                    %.5g deg\n',rmse1);
    fprintf('RMSE if half-angle scale=0.5:                 %.5g deg\n',rmse05);
    fprintf('RMSE if scale fitted:                         %.5g deg\n',rmsebest);

    % Also compare published curves with the actual MATLAB output.
    s = S{1};
    own_th = interp1(s.t,rad2deg(s.theta_b),t_eval,'linear');
    own_w = interp1(s.t,rad2deg(s.omega_b),t_eval,'linear');
    rms = @(v) sqrt(mean(v.^2));
    fprintf('\nPublished-vs-ours RMS after first 0.25 s (rough image reading):\n');
    fprintf('  Attitude, original: %.5g deg; ours/2: %.5g deg\n',...
        rms(th_pub-own_th),rms(th_pub-own_th/2));
    fprintf('  Velocity, original: %.5g deg/s; ours/2: %.5g deg/s\n',...
        rms(w_pub-own_w),rms(w_pub-own_w/2));

    figure('Color','w','Name','Xu diagnostic: published curves versus ours');
    tiledlayout(2,1,'TileSpacing','compact');
    nexttile; hold on; grid on; box on;
    plot(t_eval,w_pub,'k','LineWidth',1.4);
    plot(t_eval,own_w,'r','LineWidth',1);
    plot(t_eval,own_w/2,'b--','LineWidth',1.2);
    ylabel('Angular velocity (deg/s)');
    title('Fig. 2: published curve vs our nonzero-momenta controller');
    legend('Xu published','Ours full','Ours / 2','Location','best');
    nexttile; hold on; grid on; box on;
    plot(t_eval,th_pub,'k','LineWidth',1.4);
    plot(t_eval,own_th,'r','LineWidth',1);
    plot(t_eval,own_th/2,'b--','LineWidth',1.2);
    xlabel('Time (s)'); ylabel('Attitude (deg)');
    title('Fig. 3: published curve vs our nonzero-momenta controller');
    legend('Xu published','Ours full','Ours / 2','Location','best');

    % Compare the published paper's attitude change with integral of its
    % own published angular velocity, aligning offsets and baseline drift.
    figure('Color','w','Name','Xu paper internal rate-angle consistency');
    plot(t_eval,th_pub,'k','LineWidth',1.6); hold on; grid on; box on;
    plot(t_eval,integral_w+N*d1,'r--','LineWidth',1.3);
    plot(t_eval,0.5*integral_w+N*d05,'b:','LineWidth',1.6);
    xlabel('Time (s)'); ylabel('Attitude (degrees)');
    title('Published Fig. 3 vs integral of published Fig. 2');
    legend('Published attitude','Integrated published velocity (scale 1)',...
        'Integrated published velocity (scale 0.5)','Location','best');
    xlim([0.25 7.75]);
end

%% Helper: digitize red nonzero-momenta curve from a published PNG
function [t, y, cover] = xuDigitizeRed(filename,axisPixels,axisY,bandPixels,referenceSize)
    I = double(imread(filename));
    if size(I,3) < 3
        error('%s is not an RGB image.',filename);
    end
    [H,W,~] = size(I);
    sx = W/referenceSize(1); sy = H/referenceSize(2);
    xL = round(axisPixels(1)*sx); yTop = round(axisPixels(2)*sy);
    xR = round(axisPixels(3)*sx); yBot = round(axisPixels(4)*sy);
    yLo = round(bandPixels(1)*sy); yHi = round(bandPixels(2)*sy);
    rows = max(1,yLo):min(H,yHi);
    columns = max(1,xL):min(W,xR);
    yp = nan(numel(columns),1);
    for j = 1:numel(columns)
        col = I(rows,columns(j),:);
        r = col(:,1,1); g = col(:,1,2); b = col(:,1,3);
        good = (r>150) & (r>1.65*g) & (r>1.65*b);
        if any(good)
            selectedRows = rows(good);
            yp(j) = median(selectedRows);
        end
    end
    good = isfinite(yp);
    cover = mean(good);
    if sum(good)<5
        error('Cannot detect red plot curve in %s.',filename);
    end
    yp(~good) = interp1(columns(good),yp(good),columns(~good),'linear','extrap');
    t = (columns(:)-xL)/(xR-xL)*8;
    y = axisY(1) + (yp-yTop)/(yBot-yTop)*(axisY(2)-axisY(1));
end
