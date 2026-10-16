% Verify numerical convergence of the terminal-link inertia experiment
% by comparing post-perturbation metrics at progressively smaller timesteps.
% This tests numerical sensitivity, not physical payload capture.

clear;
clc;

% Project setup
project_root = fileparts(mfilename('fullpath'));

if ~exist(fullfile(project_root, 'config', 'paper_params.m'), 'file')
    project_root = fileparts(project_root);
end

assert(exist(fullfile(project_root, 'config', 'paper_params.m'), 'file') == 2, ...
    'Project root not found.');

addpath(genpath(project_root));

% Published step, plus two progressively smaller steps (seconds).
dts = [0.002, 0.001, 0.0005];
scales = [1, 2, 5, 10];
change_time = 4.0;

nd = numel(dts);
ns = numel(scales);
sims = cell(nd, ns);

fprintf('INERTIA PERTURBATION -- TIMESTEP REFINEMENT\n');
fprintf('Abrupt I3 change at %.1f s; results for t >= %.1f s\n\n', ...
    change_time, change_time);

for i = 1:nd
    for j = 1:ns
        cfg = paper_params();
        cfg.dt = dts(i);
        cfg.capture_time = change_time;
        cfg.payload_inertia_scale = scales(j);

        fprintf('Running dt = %.3f ms, I3 scale = %.0fx ...\n', ...
            1000*dts(i), scales(j));
        sims{i,j} = simulate_plant(cfg);

        assert(all(isfinite(sims{i,j}.theta_b(:))) && ...
               all(isfinite(sims{i,j}.delta_x(:))), ...
               'Nonfinite state detected at dt=%.4g s, scale=%.0fx', ...
               dts(i), scales(j));
    end
end

% Calculate the SAME metrics used in test_payload_capture.m.
step_ms = zeros(nd*ns,1);
inertia_scale = zeros(nd*ns,1);
peak_attitude_deg = zeros(nd*ns,1);
rms_ee_mm = zeros(nd*ns,1);
max_attitude_difference_deg = zeros(nd*ns,1);
max_tracking_difference_mm = zeros(nd*ns,1);

row = 0;
for i = 1:nd
    baseline = sims{i,1};
    t_base = baseline.t(:);

    for j = 1:ns
        s = sims{i,j};
        t = s.t(:);
        post = t >= change_time - 1e-12;
        assert(numel(t)==numel(t_base) && ...
            max(abs(t-t_base)) < 1e-10, ...
            'Time grids must match for a fixed timestep.');

        attitude_deg = rad2deg(s.theta_b(:));
        ee_mm = 1000 * vecnorm(s.delta_x, 2, 2);
        att_delta_deg = rad2deg(s.theta_b(:)-baseline.theta_b(:));
        ee_delta_mm = 1000 * vecnorm(s.delta_x-baseline.delta_x, 2, 2);

        row = row + 1;
        step_ms(row) = 1000*dts(i);
        inertia_scale(row) = scales(j);
        peak_attitude_deg(row) = max(abs(attitude_deg(post)));
        rms_ee_mm(row) = sqrt(mean(ee_mm(post).^2));
        max_attitude_difference_deg(row) = max(abs(att_delta_deg(post)));
        max_tracking_difference_mm(row) = max(ee_delta_mm(post));
    end
end

results = table(step_ms, inertia_scale, peak_attitude_deg, rms_ee_mm, ...
    max_attitude_difference_deg, max_tracking_difference_mm, ...
    'VariableNames', {'Timestep_ms','I3_Scale', ...
    'PeakAttitude_deg','RMSTracking_mm', ...
    'MaxAttitudeDifference_deg','MaxTrackingDifference_mm'});

fprintf('\n=== POST-CHANGE METRICS (same definitions as original experiment) ===\n');
disp(results);

% Quantify difference in metrics from the finest available timestep.
fine_idx = nd;
coarse_idx = 1;
mid_idx = 2;

fprintf('\n=== CHANGE IN METRICS RELATIVE TO FINEST dt=%.3f ms ===\n', ...
        1000*dts(fine_idx));
for j = 1:ns
    rcoarse = (coarse_idx-1)*ns+j;
    rmid = (mid_idx-1)*ns+j;
    rfine = (fine_idx-1)*ns+j;
    fprintf('\nI3 scale %.0fx:\n',scales(j));
    fprintf('  Peak attitude: %.3f ms = %.7f deg, %.3f ms = %.7f deg, %.3f ms = %.7f deg\n', ...
        1000*dts(coarse_idx), peak_attitude_deg(rcoarse), ...
        1000*dts(mid_idx), peak_attitude_deg(rmid), ...
        1000*dts(fine_idx), peak_attitude_deg(rfine));
    fprintf('  RMS EE error:  %.3f ms = %.7f mm,  %.3f ms = %.7f mm,  %.3f ms = %.7f mm\n', ...
        1000*dts(coarse_idx), rms_ee_mm(rcoarse), ...
        1000*dts(mid_idx), rms_ee_mm(rmid), ...
        1000*dts(fine_idx), rms_ee_mm(rfine));
    fprintf('  2-ms vs fine: |peak difference| = %.7g deg; |RMS difference| = %.7g mm\n', ...
        abs(peak_attitude_deg(rcoarse)-peak_attitude_deg(rfine)), ...
        abs(rms_ee_mm(rcoarse)-rms_ee_mm(rfine)));
end

% Compare entire attitude and EE-error trajectories at COMMON time points.
% This is a numerical convergence check, not a physical payload validation.
fine_step = sims{fine_idx,ns};
coarse_step = sims{coarse_idx,ns};
tc = coarse_step.t(:);
mask = tc >= change_time - 1e-12;
theta_fine = interp1(fine_step.t(:),fine_step.theta_b(:),tc,'linear');
dx_fine = interp1(fine_step.t(:),fine_step.delta_x,tc,'linear');
max_theta_difference = max(abs(rad2deg(coarse_step.theta_b(mask)-theta_fine(mask))));
dx_diff_mm = 1000 * vecnorm(coarse_step.delta_x(mask,:)-dx_fine(mask,:),2,2);
max_ee_vector_difference = max(dx_diff_mm);

fprintf('\n=== EXTREME CASE: 10x, 2 ms vs %.3f ms (aligned t>=4 s) ===\n',1000*dts(fine_idx));
fprintf('Max pointwise attitude difference: %.9g deg\n',max_theta_difference);
fprintf('Max pointwise EE-error VECTOR difference: %.9g mm\n',max_ee_vector_difference);

output_dir = fullfile(project_root,'results','inertia_perturbation');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
out_file = fullfile(output_dir,'inertia_timestep_convergence.csv');
writetable(results,out_file);
fprintf('\nSaved: %s\n',out_file);
fprintf('NOTE: Smaller timestep agreement supports numerical robustness only.\n');
fprintf('This remains an inertia-only step change, not a physical capture model.\n');
