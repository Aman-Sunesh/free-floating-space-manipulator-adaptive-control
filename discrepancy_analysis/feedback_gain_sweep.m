% Compare simulated spacecraft angular velocity and attitude against
% digitized red curves from Xu's Figs. 2 and 3.
% Sweep the attitude-feedback gain to investigate the reproduction mismatch.
% Screenshot measurements are approximate, not original numerical data.
% Evaluate the interval 0.25-7.75 s to exclude endpoint transients.

addpath(genpath(pwd));
assert(isfile('actual_fig2.png') && isfile('actual_fig3.png'), ...
       'Put actual_fig2.png and actual_fig3.png in the current folder first.');

cfg = paper_params();
t_eval = (0.25:0.005:7.75).';
[tw,wp] = digitize_published_red('actual_fig2.png', [136 32 1061 763], ...
    [1 -3.5], [141 250], [1118 855]);
[ta,ap] = digitize_published_red('actual_fig3.png', [199 23 1361 939], ...
    [0.015 -0.030], [170 750], [1403 1032]);
wp = movmean(interp1(tw,wp,t_eval,'linear'),11);
ap = movmean(interp1(ta,ap,t_eval,'linear'),11);

mult = [0.5 0.75 1 1.25 1.5 1.75 2 2.5 3];
rmse_w = zeros(size(mult));
rmse_a = zeros(size(mult));
fit_w = zeros(size(mult));
fit_a = zeros(size(mult));
store = cell(size(mult));

fprintf('\nGAIN SWEEP (published red curve vs simulated physical outputs)\n');
fprintf('Gain x | lambda_b | RMSE omega (deg/s) | RMSE theta (deg) | LS amplitude factors omega,theta\n');
for k = 1:numel(mult)
    c = cfg;

    % Scale attitude feedback in both the controller and adaptation laws.
    c.lambda_b = cfg.lambda_b*mult(k);

    s = simulate_plant(c);
    w = interp1(s.t,rad2deg(s.omega_b),t_eval,'linear');
    a = interp1(s.t,rad2deg(s.theta_b),t_eval,'linear');
    store{k} = struct('omega',w,'theta',a);
    rmse_w(k) = sqrt(mean((wp-w).^2));
    rmse_a(k) = sqrt(mean((ap-a).^2));
    fit_w(k) = (w'*wp)/(w'*w);
    fit_a(k) = (a'*ap)/(a'*a);
    fprintf('%5.2f  | %8.2f | %18.7f | %16.7f | %7.3f, %7.3f\n', ...
        mult(k),c.lambda_b,rmse_w(k),rmse_a(k),fit_w(k),fit_a(k));
end

% Report winners separately; do not combine unlike physical units arbitrarily.
[~,kw] = min(rmse_w);
[~,ka] = min(rmse_a);
score = (rmse_w/rmse_w(mult==1)).^2 + (rmse_a/rmse_a(mult==1)).^2;
[~,kc] = min(score);
fprintf('\nBest omega match:       gain x%.2f\n',mult(kw));
fprintf('Best attitude match:    gain x%.2f\n',mult(ka));
fprintf('Best combined match:    gain x%.2f (equal relative weighting)\n',mult(kc));
fprintf('Baseline gain x1 has omega RMSE %.7g and attitude RMSE %.7g\n', ...
    rmse_w(mult==1),rmse_a(mult==1));

figure('Color','w','Name','Physical feedback gain vs published figure curves');
tiledlayout(2,1,'TileSpacing','compact');
nexttile; hold on; grid on; box on;
plot(t_eval,wp,'k','LineWidth',1.8);
plot(t_eval,store{find(mult==1,1)}.omega,'r:','LineWidth',1.3);
plot(t_eval,store{kc}.omega,'b--','LineWidth',1.4);
ylabel('Angular velocity (deg/s)');
legend('Xu published', 'Our original gain',sprintf('Best gain x%.2f',mult(kc)), 'Location','best');
title('Published angular velocity vs actual simulated output');
nexttile; hold on; grid on; box on;
plot(t_eval,ap,'k','LineWidth',1.8);
plot(t_eval,store{find(mult==1,1)}.theta,'r:','LineWidth',1.3);
plot(t_eval,store{kc}.theta,'b--','LineWidth',1.4);
ylabel('Attitude (deg)'); xlabel('Time (s)');
legend('Xu published','Our original gain',sprintf('Best gain x%.2f',mult(kc)), 'Location','best');
title('Published attitude vs actual simulated output');

figure('Color','w','Name','Gain sweep quantitative errors');
tiledlayout(2,1,'TileSpacing','compact');
nexttile; plot(mult,rmse_w,'o-','LineWidth',1.5); grid on; box on;
ylabel('Angular velocity RMSE (deg/s)');
nexttile; plot(mult,rmse_a,'o-','LineWidth',1.5); grid on; box on;
xlabel('Controller attitude gain multiplier'); ylabel('Attitude RMSE (deg)');

function [t,y] = digitize_published_red(filename,ax,ay,bands,sz)
    I = double(imread(filename));
    [H,W,~] = size(I); sx=W/sz(1); sy=H/sz(2);
    xL=round(ax(1)*sx); xR=round(ax(3)*sx);
    yTop=round(ax(2)*sy); yBot=round(ax(4)*sy);
    cols=max(1,xL):min(W,xR);
    rows=max(1,round(bands(1)*sy)):min(H,round(bands(2)*sy));
    yp=nan(numel(cols),1);
    for j=1:numel(cols)
        rgb=I(rows,cols(j),:);
        r=rgb(:,1,1); g=rgb(:,1,2); b=rgb(:,1,3);
        good=(r>150) & (r>1.65*g) & (r>1.65*b);
        if any(good), yp(j)=median(rows(good)); end
    end
    ok=isfinite(yp);
    assert(mean(ok)>0.8,'Published red curve extraction failed for %s',filename);
    yp(~ok)=interp1(cols(ok),yp(ok),cols(~ok),'linear','extrap');
    t=(cols(:)-xL)/(xR-xL)*8;
    y=ay(1)+(yp-yTop)/(yBot-yTop)*(ay(2)-ay(1));
end
