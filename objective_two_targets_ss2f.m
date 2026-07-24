function loss = objective_two_targets_ss2f(u, p, target_exit, target_rel_size_H)

%OBJECTIVE_TWO_TARGETS_SS2F
%Objective for two-type steady-state calibration
%Lambda and gamma are calibrated from data

%Calibrated:
%u(1) = log(x)
%u(2) = log(A_H)

%Targets:
%x   -> aggregate exit rate
%A_H -> relative BPW/output per firm of high energy-intensive firms

x   = exp(u(1));
A_H = exp(u(2));

if x <= 0 || A_H <= 0 || ~isfinite(x) || ~isfinite(A_H)
    loss = 1e12;
    return;
end

try
    out = SS_moments_twotype_AH_x(A_H, x, p);

    model_rel_size_H = (out.Y_type(2) / out.M_type(2)) / ...
                       (out.Y_type(1) / out.M_type(1));

    err_exit = (out.exit_rate_agg - target_exit) / max(target_exit, 1e-8);
    err_size = (model_rel_size_H - target_rel_size_H) / max(target_rel_size_H, 1e-8);

    loss = err_exit^2 + err_size^2;

    if ~isfinite(loss)
        loss = 1e12;
    end

catch ME
    warning('objective_two_targets_ss2f failed: %s', ME.message);
    loss = 1e12;
end

end
