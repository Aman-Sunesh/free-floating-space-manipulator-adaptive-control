% Final verification of Xu's adaptive controller reproduction.
% Compare 1x and 2x attitude gains for both controller variants,
% examine the initial 0.25 s, and verify the adaptation updates.
% Published curves are approximate measurements from screenshots.

clc; close all;
addpath(genpath(pwd));
assert(exist('paper_params','file')==2, 'Start MATLAB in the reproduction project root.');
assert(isfile('actual_fig2.png') && isfile('actual_fig3.png'), ...
    'Copy actual_fig2.png and actual_fig3.png to this folder.');

cfg = paper_params();
gains = [1 2];
S = cell(2,2); % rows: gain 1x, 2x; cols: nonzero, zero mom.
for ig=1:2
    c = cfg; c.lambda_b = cfg.lambda_b*gains(ig);
    fprintf('\nRunning actual controller at gain x%g ...\n',gains(ig));
    S{ig,1} = simulate_plant(c);
    S{ig,2} = simulate_zero_momentum_plant(c);
end

%% 1. Published comparison: red and blue spacecraft ATTITUDE
% Calibrations match the original uploaded actual_fig2/3 PNG files.
t_eval=(0.25:0.005:7.75).';
[taR, thPubR, covR]=readXuTrace('actual_fig3.png', [199 23 1361 939], ...
    [0.015 -0.030], [170 750], [1403 1032], 'red');
[taB, thPubB, covB]=readXuTrace('actual_fig3.png', [199 23 1361 939], ...
    [0.015 -0.030], [125 920], [1403 1032], 'blue');
[twR, wPubR, covWR]=readXuTrace('actual_fig2.png', [136 32 1061 763], ...
    [1 -3.5], [141 250], [1118 855], 'red');
% A blue Fig.2 extraction is supplemental only: dashed image trace is noisy.
[twB, wPubB, covWB]=readXuTrace('actual_fig2.png', [136 32 1061 763], ...
    [1 -3.5], [150 245], [1118 855], 'blue');

thPubR = movmean(interp1(taR,thPubR,t_eval,'linear'),11);
thPubB = movmean(interp1(taB,thPubB,t_eval,'linear'),11);
wPubR = movmean(interp1(twR,wPubR,t_eval,'linear'),11);
wPubB = movmean(interp1(twB,wPubB,t_eval,'linear'),11);
rmsf = @(v) sqrt(mean(v.^2));

fprintf('\n==================== TEST 1: ZERO-MOMENTA GAIN ====================\n');
fprintf('Published image coverage: red Fig3 %.1f%%, blue Fig3 %.1f%%, red Fig2 %.1f%%, blue Fig2 %.1f%%\n', ...
    100*covR,100*covB,100*covWR,100*covWB);
ERR = zeros(2,4); % red theta, red omega, blue theta, blue omega
for ig=1:2
    a = S{ig,1}; b = S{ig,2};
    an=interp1(a.t,rad2deg(a.theta_b),t_eval,'linear');
    wn=interp1(a.t,rad2deg(a.omega_b),t_eval,'linear');
    az=interp1(b.t,rad2deg(b.theta_b),t_eval,'linear');
    wz=interp1(b.t,rad2deg(b.omega_b),t_eval,'linear');
    ERR(ig,:)=[rmsf(thPubR-an),rmsf(wPubR-wn), ...
                rmsf(thPubB-az),rmsf(wPubB-wz)];
    fprintf('Gain x%g: NONZERO theta RMSE %.7g deg, omega %.7g deg/s; ZERO theta %.7g deg, omega %.7g deg/s\n', ...
        gains(ig),ERR(ig,1),ERR(ig,2),ERR(ig,3),ERR(ig,4));
end
fprintf('Improvement ratio x1/x2: nonzero attitude %.2fx, omega %.2fx; zero attitude %.2fx, omega %.2fx\n', ...
    ERR(1,1)/ERR(2,1),ERR(1,2)/ERR(2,2),ERR(1,3)/ERR(2,3),ERR(1,4)/ERR(2,4));
if ERR(2,3)<ERR(1,3)
    fprintf('PASS: 2x improves ZERO-MOMENTA published attitude agreement.\n');
else
    fprintf('NOTE: 2x does NOT improve zero-momenta attitude. Do not claim it does.\n');
end
fprintf('NOTE: blue Fig2 is dashed/noisy; treat its numeric RMSE as exploratory.\n');

figure('Color','w','Name','Final test: original vs double gain, BOTH controllers');
tiledlayout(2,2,'TileSpacing','compact');
for j=1:2
    if j==1, tp=thPubR; wp=wPubR; name='Nonzero-momenta';
    else, tp=thPubB; wp=wPubB; name='Zero-momenta'; end
    nexttile; hold on; grid on;
    plot(t_eval,tp,'k','LineWidth',1.7);
    for ig=1:2
        a=S{ig,j};
        plot(t_eval,interp1(a.t,rad2deg(a.theta_b),t_eval), ...
            pickstyle(ig),'LineWidth',1.2);
    end
    title([name,' attitude']); ylabel('deg'); xlabel('Time (s)');
    legend('Published','x1','x2','Location','best');
    nexttile; hold on; grid on;
    plot(t_eval,wp,'k','LineWidth',1.7);
    for ig=1:2
        a=S{ig,j};
        plot(t_eval,interp1(a.t,rad2deg(a.omega_b),t_eval), ...
            pickstyle(ig),'LineWidth',1.2);
    end
    title([name,' angular velocity']); ylabel('deg/s'); xlabel('Time (s)');
    legend('Published','x1','x2','Location','best');
end
exportgraphics(gcf,'Xu_FINAL_gain_both_controllers.png','Resolution',200);

%% 2. Early transient and momentum: no unreliable pixel fitting at t<0.25
fprintf('\n==================== TEST 2: FIRST 0.25 SECONDS ================\n');
preOmega=rad2deg(cfg.qb_dot_0);
fprintf('Published prescribed PRE-CONTROL initial omega: %.9f deg/s\n',preOmega);
for j=1:2
    if j==1, name='Nonzero-momenta'; else, name='Zero-momenta'; end
    fprintf('\n%s:\n',name);
    for ig=1:2
        a=S{ig,j};
        ix=a.t<=0.25;
        dt=diff(a.t);
        rateResidual=diff(a.theta_b)./dt-a.omega_b(1:end-1);
        [Hb,Hbm]=momentum_matrices(a.phi(1,:).',cfg.l,cfg.r,cfg.m,cfg.I);
        wMomentum=(a.p0_true-Hbm*a.phi_dot(1,:).')/Hb;
        maxW=max(abs(rad2deg(a.omega_b(ix))));
        maxTh=max(abs(rad2deg(a.theta_b(ix))));
        fprintf('  x%g first LOGGED omega %.9f deg/s; predicted %.9f; early max|omega| %.7g; early max|theta| %.7g\n', ...
            gains(ig),rad2deg(a.omega_b(1)),rad2deg(wMomentum),maxW,maxTh);
        fprintf('       max Euler |dtheta/dt - omega| = %.3e deg/s\n', ...
            max(abs(rad2deg(rateResidual(1:sum(ix)-1)))));
        assert(max(abs(rateResidual))<1e-8,'Failed dtheta/dt=omega check.');
        assert(abs(a.omega_b(1)-wMomentum)<1e-10,'Failed initial momentum check.');
    end
    fprintf('  gain x1 vs x2 first logged omega difference: %.3e deg/s\n', ...
        abs(rad2deg(S{1,j}.omega_b(1)-S{2,j}.omega_b(1))));
end
fprintf('NOTE: Fig2 has a red annotation and vertical spike at t~0, so an image-based\n');
fprintf('      RMSE of the first 0.25 s cannot be trusted without raw published data.\n');
figure('Color','w','Name','First quarter-second and controller-gain comparison');
tiledlayout(2,1,'TileSpacing','compact');
nexttile; hold on; grid on;
for j=1:2
  for ig=1:2
    s=S{ig,j}; i=s.t<=0.25;
    plot(s.t(i),rad2deg(s.omega_b(i)),'LineWidth',1.2, ...
        'LineStyle',pickstyle2(ig),'Color',pickcolor(j));
  end
end
xlabel('Time (s)'); ylabel('omega (deg/s)');
title('First 0.25 seconds (actual simulation, not extracted image)');
legend('nonzero x1','nonzero x2','zero x1','zero x2','Location','best');
nexttile; hold on; grid on;
for j=1:2
  for ig=1:2
    s=S{ig,j}; i=s.t<=0.25;
    plot(s.t(i),rad2deg(s.theta_b(i)),'LineWidth',1.2, ...
        'LineStyle',pickstyle2(ig),'Color',pickcolor(j));
  end
end
xlabel('Time (s)'); ylabel('attitude (deg)');
legend('nonzero x1','nonzero x2','zero x1','zero x2','Location','best');
exportgraphics(gcf,'Xu_FINAL_first_025s.png','Resolution',200);

%% 3. Verify adaptation laws and Euler parameter updates
% Xu et al. (2016), Eqs. (12), (14), (22), and (23).
% Reconstruct the updates without calling parameter_update.m.
fprintf('\n==================== TEST 3: ADAPTATION LAWS ====================\n');
allOK=true;
for ig=1:2
  c=cfg; c.lambda_b=cfg.lambda_b*gains(ig);
  for j=1:2
    s=S{ig,j};
    idx=unique(round(linspace(1,numel(s.t)-1,120)));
    maxDynamic=0; maxKinematic=0; maxY1=0; maxY2=0;
    maxMomentum=0; maxPhiLaw=0;
    for k=idx(:).'
      h=s.t(k+1)-s.t(k);
      theta=s.theta_b(k); phi=s.phi(k,:).'; x=s.x(k,:).';
      omega=s.omega_b(k); phiDot=s.phi_dot(k,:).';
      thetaError=atan2(sin(theta-c.qbd),cos(theta-c.qbd));
      epsBV=sin(thetaError/2);  % checked quaternion convention
      [xd,xdDot]=desired_trajectory(s.t(k));
      xDot=(s.x(k+1,:).'-x)/h;
      dx=x-xd; dxDot=xDot-xdDot;
      if j==1
        adHat=s.ad_bar_hat(k,1:10).'; akHat=s.ak_bar_hat(k,1:4).';
        pHat=s.ad_bar_hat(k,11); vHat=s.ak_bar_hat(k,5:6).';
        [YdBar,Yd]=dynamic_regressor(phi,omega,phiDot);
        [YkBar,Yk]=kinematic_regressor(theta,phi,omega,phiDot);
        gd=c.Gamma_d; gk=c.Gamma_k;
        adObserved=((s.ad_bar_hat(k+1,:)-s.ad_bar_hat(k,:))/h).';
        akObserved=((s.ak_bar_hat(k+1,:)-s.ak_bar_hat(k,:))/h).';
      else
        adHat=s.ad_hat(k,:).'; akHat=s.ak_hat(k,:).';
        pHat=0; vHat=zeros(2,1);
        [~,Yd]=dynamic_regressor(phi,omega,phiDot);
        [~,Yk]=kinematic_regressor(theta,phi,omega,phiDot);
        YdBar=Yd; YkBar=Yk;
        gd=c.Gamma_d(1:10,1:10); gk=c.Gamma_k(1:4,1:4);
        adObserved=((s.ad_hat(k+1,:)-s.ad_hat(k,:))/h).';
        akObserved=((s.ak_hat(k+1,:)-s.ak_hat(k,:))/h).';
      end
      [~,YHb]=dynamic_regressor(phi,1,zeros(3,1));
      HbHat=YHb*adHat;
      HbmHat=zeros(1,3); JmHat=zeros(2,3);
      for q=1:3
        e=zeros(3,1); e(q)=1;
        [~,Yhbm]=dynamic_regressor(phi,0,e);
        HbmHat(q)=Yhbm*adHat;
        [~,Yjm]=kinematic_regressor(theta,phi,0,e);
        JmHat(:,q)=Yjm*akHat;
      end
      [~,Yjb]=kinematic_regressor(theta,phi,1,zeros(3,1));
      JbHat=Yjb*akHat;
      Th=eye(3)-pinv(HbmHat)*HbmHat;
      phiAtt=pinv(HbmHat)*(pHat+HbHat*c.lambda_b*epsBV);
      remainder=-vHat+xdDot-c.Lambda_x*dx-JmHat*phiAtt;
      phiCmd=Th*(pinv(JmHat*Th)*remainder)+phiAtt;
      % Ideal execution: phiDot = phiCmd.
      y1=HbHat*(omega+c.lambda_b*epsBV)+HbmHat*(phiDot-phiCmd);
      y2=JbHat*omega+JmHat*(phiDot-phiCmd)-(dxDot+c.Lambda_x*dx);
      adExpected=-gd*(YdBar.')*y1;
      akExpected=-gk*(YkBar.')*y2;
      maxDynamic=max(maxDynamic,norm(adObserved-adExpected,inf));
      maxKinematic=max(maxKinematic,norm(akObserved-akExpected,inf));
      maxY1=max(maxY1,abs(y1-s.y1(k)));
      maxY2=max(maxY2,norm(y2-s.y2(k,:).',inf));
      maxPhiLaw=max(maxPhiLaw,norm(phiDot-phiCmd,inf));
      [Hb,Hbm]=momentum_matrices(phi,c.l,c.r,c.m,c.I);
      maxMomentum=max(maxMomentum,abs(Hb*omega+Hbm*phiDot-s.p0_true));
    end
    name={'nonzero','zero'};
    thisOK=(maxDynamic<1e-7)&&(maxKinematic<1e-7)&& ...
           (maxY1<1e-7)&&(maxY2<1e-7)&& ...
           (maxPhiLaw<1e-7)&&(maxMomentum<1e-8);
    allOK=allOK && thisOK;
    if thisOK, status='PASS'; else, status='CHECK'; end
    fprintf('%s: %s controller, gain x%g, %d sampled instants\n', ...
        status,name{j},gains(ig),numel(idx));
    fprintf('  max residual: dynamic update %.3e, kinematic update %.3e\n',maxDynamic,maxKinematic);
    fprintf('  max residual: y1 %.3e, y2 %.3e, commanded joint velocity %.3e, momentum %.3e\n', ...
        maxY1,maxY2,maxPhiLaw,maxMomentum);
  end
end
if allOK
 fprintf('\nPASS: both controllers at both gains obey the implemented Xu adaptation laws.\n');
else
 warning('At least one equation check failed. Investigate before concluding.');
end
fprintf('Saved: Xu_FINAL_gain_both_controllers.png and Xu_FINAL_first_025s.png\n');

%% Helpers - used only for plots and image digitization
function style=pickstyle(i)
  if i==1,style='r:';else,style='b--';end
end
function style=pickstyle2(i)
  if i==1,style='-';else,style='--';end
end
function color=pickcolor(j)
  if j==1,color=[0.9,0,0];else,color=[0,0.25,0.8];end
end
function [t,y,coverage]=readXuTrace(filename,ax,ay,bands,refSize,which)
  A=double(imread(filename));
  [H,W,~]=size(A); sx=W/refSize(1); sy=H/refSize(2);
  xL=round(ax(1)*sx); xR=round(ax(3)*sx);
  yT=round(ax(2)*sy); yB=round(ax(4)*sy);
  col=max(1,xL):min(W,xR);
  row=max(1,round(bands(1)*sy)):min(H,round(bands(2)*sy));
  yp=nan(numel(col),1);
  for ii=1:numel(col)
    rgb=A(row,col(ii),:); R=rgb(:,1,1);G=rgb(:,1,2);B=rgb(:,1,3);
    if strcmp(which,'red')
      ok=(R>150)&(R>1.65*G)&(R>1.65*B);
    else
      ok=(B>R+22)&(B>G+18)&(B>100)&(R<180);
    end
    if any(ok),yp(ii)=median(row(ok));end
  end
  valid=isfinite(yp); coverage=mean(valid);
  assert(coverage>0.4,'Cannot reliably digitize %s from %s.',which,filename);
  yp(~valid)=interp1(col(valid),yp(valid),col(~valid),'linear','extrap');
  t=(col(:)-xL)/(xR-xL)*8;
  y=ay(1)+(yp-yT)/(yB-yT)*(ay(2)-ay(1));
end
