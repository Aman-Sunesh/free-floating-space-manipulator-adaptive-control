clear;
clc;
close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

% Test robustness to imperfect joint-velocity execution.
% alpha = 1: ideal execution; alpha < 1: slower; alpha > 1: faster.
alphas = [1.2 1.1 1.05 1.0 0.98 0.95 0.9 0.8];

sims = cell(size(alphas));

fprintf('========================================\n');
fprintf('JOINT-VELOCITY SERVO MISMATCH TEST\n');
fprintf('========================================\n\n');

for k = 1:length(alphas)

    cfg = paper_params();

    % Scale the executed joint velocity while retaining the same controller.
    cfg.servo_alpha = alphas(k);

    fprintf('alpha = %.2f\n', alphas(k));

    sims{k} = simulate_plant(cfg);

end


%% Plot spacecraft attitude
figure('Color','w');
hold on;
box on;

for k = 1:length(alphas)

    plot( ...
        sims{k}.t, ...
        rad2deg(sims{k}.theta_b), ...
        'LineWidth', 1.2);

end

xlabel('Time (s)');
ylabel('Spacecraft attitude (deg)');
legend(compose('\\alpha = %.2f', alphas), ...
    'Location','best');

title('Effect of Joint-Velocity Servo Mismatch');


%% Plot end-effector tracking error
figure('Color','w');
hold on;
box on;

for k = 1:length(alphas)

    tracking_error = vecnorm( ...
        sims{k}.delta_x, 2, 2);

    plot( ...
        sims{k}.t, ...
        tracking_error, ...
        'LineWidth', 1.2);

end

xlabel('Time (s)');
ylabel('EE tracking-error norm (m)');
legend(compose('\\alpha = %.2f', alphas), ...
    'Location','best');

title('Effect on End-Effector Tracking');


%% Print summary
fprintf('\n');
fprintf('alpha    max attitude (deg)    final EE error (m)\n');
fprintf('-------------------------------------------------\n');

for k = 1:length(alphas)

    max_att = max(abs( ...
        rad2deg(sims{k}.theta_b)));

    final_err = norm( ...
        sims{k}.delta_x(end,:));

    fprintf('%5.2f       %10.4f          %.6e\n', ...
        alphas(k), max_att, final_err);

end