% Two-type steady states under Tauchen and Rouwenhorst discretizations.
% The calibration matches aggregate energy use, exit, and relative firm size.
clear; clc;

code_dir = fileparts(mfilename('fullpath'));
root = fileparts(code_dir);
addpath(code_dir, fullfile(code_dir, 'helpers'));
baselines_dir = fullfile(root, 'data', 'baselines');
calibration_dir = fullfile(root, 'data', 'calibration');

common_target_source = 'tauchen'; % 'tauchen' or 'rouwenhorst'

switch lower(common_target_source)
    case 'tauchen'
        load(fullfile(baselines_dir, ...
            'German_Baseline_OneType_TAUCHEN.mat'),'e_agg');
    case 'rouwenhorst'
        load(fullfile(baselines_dir, ...
            'German_Baseline_OneType_ROUWENHORST.mat'),'e_agg');
    otherwise
        error('Unknown common_target_source.');
end

common_target_e_agg = e_agg;

methods = {'tauchen','rouwenhorst'};
gamma_rule = 'median_split'; % 'median_split' baseline, 'tercile_split', 'quartile_split'

%% =========================================================
% 1. READ GAMMA AND LAMBDA FROM GERMAN DATA
% ==========================================================
% gamma_L and gamma_H come from energy intensity:
%   gamma_s = Energieverbrauch_s / Bruttoproduktionswert_s
%
% lambda_L and lambda_H come from firm counts:
%   lambda_s = Unternehmen_s / Unternehmen_total
%
% Both are computed using the same sector classification rule.

gamma_file = fullfile(calibration_dir, 'gamma_classification_rules.csv');

if exist(gamma_file, 'file') == 2

    Gtab = readtable(gamma_file, 'TextType', 'string');

    ridx = find(strcmpi(strtrim(string(Gtab.rule)), gamma_rule), 1);

    if isempty(ridx)
        error('gamma_rule %s not found in %s', gamma_rule, gamma_file);
    end

    gamma_L_data = Gtab.gamma_L(ridx);
    gamma_H_data = Gtab.gamma_H(ridx);

    % IMPORTANT:
    % lambda_H must be the firm-count share, not share_high_bpw.
    if any(strcmpi(Gtab.Properties.VariableNames,'lambda_H'))
        lambda_H_data = Gtab.lambda_H(ridx);
    else
        error(['Column lambda_H not found in %s.\n' ...
               'Re-run calibrate_gamma_from_sector_data.m with the Unternehmen column included.'], gamma_file);
    end

    if any(strcmpi(Gtab.Properties.VariableNames,'lambda_L'))
        lambda_L_data = Gtab.lambda_L(ridx);
    else
        lambda_L_data = 1 - lambda_H_data;
    end

    if any(strcmpi(Gtab.Properties.VariableNames,'rel_bpw_per_firm_H'))
        target_rel_size_H = Gtab.rel_bpw_per_firm_H(ridx);
    else
        error(['Column rel_bpw_per_firm_H not found in %s.\n' ...
            'Re-run calibrate_gamma_from_sector_data.m after adding the BPW-per-firm calculation.'], gamma_file);
    end

    fprintf('Using classification rule %s from %s\n', gamma_rule, gamma_file);
    fprintf('gamma_L  = %.8f\n', gamma_L_data);
    fprintf('gamma_H  = %.8f\n', gamma_H_data);
    fprintf('lambda_L = %.8f\n', lambda_L_data);
    fprintf('lambda_H = %.8f\n', lambda_H_data);

else

    error(['%s not found.\n' ...
           'First run calibrate_gamma_from_sector_data.m to create gamma_classification_rules.csv.'], gamma_file);

end

%% =========================================================
% 2. LOOP OVER DISCRETIZATION METHODS
% ==========================================================
for m = 1:numel(methods)

    method = methods{m};

    fprintf('\n====================================================\n');
    fprintf('Running two-type baseline with %s path\n', upper(method));
    fprintf('====================================================\n');

    target_e_agg = common_target_e_agg;

    death_rate_2018 = 0.033557128;
    death_rate_2019 = 0.033450969;

    target_exit = 0.5 * (death_rate_2018 + death_rate_2019);

    %% =====================================================
    % 3. FIXED PARAMETERS
    % ======================================================
    p = struct();

    p.beta        = 0.95^(1/12);
    p.rho_theta   = 0.9949;
    p.sigma_theta = 0.0564;

    p.alpha    = 0.22;
    p.L_supply = 1;

    p.W  = 1;
    p.Pe = 1;

    p.n_theta = 101;

    % Data-calibrated energy intensities
    p.gamma_type = [gamma_L_data, gamma_H_data];

    % Data-calibrated firm-count type shares
    lambda_H0 = lambda_H_data;
    lambda_L0 = lambda_L_data;

    p.lambda_type = [lambda_L0, lambda_H0];

    % Store data moments for reporting
    p.lambda_L_data = lambda_L_data;
    p.lambda_H_data = lambda_H_data;

    %% =====================================================
    % 4. PRODUCTIVITY GRID AND TRANSITION MATRIX
    % ======================================================
    switch lower(method)
        case 'tauchen'
            [z, p.Q_incumbent] = tauchen(p.n_theta, 0, p.rho_theta, p.sigma_theta, 4);

        case 'rouwenhorst'
            [z, p.Q_incumbent] = rouwenhorst_ar1(p.n_theta, 0, p.rho_theta, p.sigma_theta);
    end

    % z is log productivity, theta is level productivity
    p.theta_grid = exp(z);

    % Entrant productivity distribution: stationary distribution of Q_incumbent
    [V,D] = eig(p.Q_incumbent');
    [~,idx] = min(abs(diag(D) - 1));

    pi_stat = real(V(:,idx));
    pi_stat = pi_stat / sum(pi_stat);

    p.Q_entrant = pi_stat;

    clear V D idx pi_stat z

   %% =====================================================
% 5. CALIBRATE x AND A_H ONLY
% ======================================================
% gamma_L and gamma_H are fixed from sector-level energy intensity data.
% lambda_L and lambda_H are fixed from firm-count shares.
%
% Internally calibrated parameters:
%   x   -> aggregate exit rate
%   A_H -> relative BPW/output per firm of high-energy firms
%
% The low-energy productivity shifter is normalized to A_L = 1.
% A_H is the high-type productivity/scale shifter. Since high-energy
% firms differ from low-energy firms not only in energy intensity but also
% in average firm size, A_H allows the model to match the empirical
% relative BPW per firm of high-energy sectors.
%
% Aggregate energy use is not directly targeted in the objective.
% It is used as a validation moment after the steady state is solved.

    obj = @(u) objective_two_targets_ss2f(u, p, target_exit, target_rel_size_H);

    opts = optimset('Display','iter', ...
                    'TolX',1e-8, ...
                    'TolFun',1e-10, ...
                    'MaxIter',500, ...
                    'MaxFunEvals',2500);

    % u(1) = log(x)
    % u(2) = log(A_H)
    starts = [ ...
        log(0.500), log(1.00); ...
        log(0.450), log(0.95); ...
        log(0.550), log(1.05); ...
        log(0.500), log(1.15); ...
        log(0.600), log(0.90); ...
        log(0.400), log(1.25)];

    best_loss = inf;
    ustar = starts(1,:)';

    for si = 1:size(starts,1)

        [u_try, f_try] = fminsearch(obj, starts(si,:)', opts);

        if f_try < best_loss
            best_loss = f_try;
            ustar = u_try;
        end

    end

    x_cal   = exp(ustar(1));
    A_H_cal = exp(ustar(2));

    %% =====================================================
    % 6. SOLVE STEADY STATE
    % ======================================================
    out = SS_moments_twotype_AH_x(A_H_cal, x_cal, p);
    model_rel_size_H = (out.Y_type(2) / out.M_type(2)) / ...
        (out.Y_type(1) / out.M_type(1));

    err_size = model_rel_size_H - target_rel_size_H;

    err_e = out.e_agg - target_e_agg;
    err_x = out.exit_rate_agg - target_exit;

    %% =====================================================
    % 7. STORE CALIBRATED PARAMETERS
    % ======================================================
    gamma_L    = p.gamma_type(1);
    gamma_H    = p.gamma_type(2);
    gamma_type = p.gamma_type;

    lambda_type = p.lambda_type;
    lambda_L    = lambda_type(1);
    lambda_H    = lambda_type(2);

    A_L    = 1.0;
    A_H    = A_H_cal;
    A_type = [A_L, A_H];

    alpha       = p.alpha;
    beta        = p.beta;
    rho_theta   = p.rho_theta;
    sigma_theta = p.sigma_theta;

    L_supply = p.L_supply;
    W        = p.W;
    Pe       = p.Pe;

    n_theta     = p.n_theta;
    theta_grid  = p.theta_grid;
    Q_incumbent = p.Q_incumbent;
    Q_entrant   = p.Q_entrant;

    x  = x_cal;
    xE = out.xE_implied;

    V_type      = out.V_type;
    exit_policy = out.exit_policy;

    mu         = out.mu;
    mu_operate = out.mu_operate;

    share_type     = out.share_type;
    exit_rate_type = out.exit_rate_type;
    exit_rate_agg  = out.exit_rate_agg;

    M_total          = out.M_total;
    M_type           = out.M_type;
    M_entrants_total = out.M_entrants_total;
    M_entrants_type  = out.M_entrants_type;

    y_policy  = out.y_policy;
    n_policy  = out.n_policy;
    e_policy  = out.e_policy;
    pi_policy = out.pi_policy;

    y_agg  = out.y_agg;
    n_agg  = out.n_agg;
    e_agg  = out.e_agg;
    pi_agg = NaN;

    Y_type = out.Y_type;
    N_type = out.N_type;
    E_type = out.E_type;

    output_share_type = out.output_share_type;
    labor_share_type  = out.labor_share_type;
    energy_share_type = out.energy_share_type;

    %% =====================================================
    % 8. SAVE BASELINE
    % ======================================================
    fname = fullfile(baselines_dir, ...
        sprintf('German_Baseline_TwoTypes_%s.mat', upper(method)));

    save(fname, ...
        'method','gamma_rule','common_target_source', ...
        'target_e_agg','target_exit', ...
        'beta','rho_theta','sigma_theta','alpha', ...
        'x','xE', ...
        'gamma_L','gamma_H','gamma_type', ...
        'lambda_L','lambda_H','lambda_type', ...
        'A_L','A_H','A_type', ...
        'W','Pe','L_supply', ...
        'n_theta','theta_grid','Q_incumbent','Q_entrant', ...
        'V_type','exit_policy','mu','mu_operate', ...
        'share_type','exit_rate_type','exit_rate_agg', ...
        'M_total','M_type','M_entrants_total','M_entrants_type', ...
        'y_policy','n_policy','e_policy','pi_policy', ...
        'y_agg','n_agg','e_agg','pi_agg', ...
        'Y_type','N_type','E_type', ...
        'output_share_type','labor_share_type','energy_share_type', 'target_rel_size_H','model_rel_size_H');

    %% =====================================================
    % 9. PRINT RESULTS
    % ======================================================
    fprintf('\nSaved: %s\n', fname);

    fprintf('\nTwo-type calibration results\n');
    fprintf('----------------------------\n');
    fprintf('x calibrated        = %.8f\n', x_cal);

    fprintf('Using classification rule %s from %s\n', gamma_rule, gamma_file);
    fprintf('gamma_L  = %.8f\n', gamma_L_data);
    fprintf('gamma_H  = %.8f\n', gamma_H_data);
    fprintf('lambda_L = %.8f\n', lambda_L_data);
    fprintf('lambda_H = %.8f\n', lambda_H_data);
    fprintf('target relative size H = %.8f\n', target_rel_size_H);

    fprintf('\nTargets and model moments\n');
    fprintf('-------------------------\n');

    fprintf('Target exit              = %.8f\n', target_exit);
    fprintf('Model exit               = %.8f\n', out.exit_rate_agg);
    fprintf('Exit error               = %.3e\n', err_x);

    fprintf('Target rel size H        = %.8f\n', target_rel_size_H);
    fprintf('Model rel size H         = %.8f\n', model_rel_size_H);
    fprintf('Rel size error           = %.3e\n', err_size);

    fprintf('\nValidation moment\n');
    fprintf('-----------------\n');
    fprintf('Target e_agg             = %.8f\n', target_e_agg);
    fprintf('Model e_agg              = %.8f\n', out.e_agg);
    fprintf('Energy validation error  = %.3e\n', err_e);

    fprintf('\nType moments\n');
    fprintf('------------\n');
    fprintf('Active share L      = %.8f\n', out.share_type(1));
    fprintf('Active share H      = %.8f\n', out.share_type(2));
    fprintf('Exit rate L         = %.8f\n', out.exit_rate_type(1));
    fprintf('Exit rate H         = %.8f\n', out.exit_rate_type(2));
    fprintf('Energy share L      = %.8f\n', out.energy_share_type(1));
    fprintf('Energy share H      = %.8f\n', out.energy_share_type(2));
    fprintf('A_H calibrated      = %.8f\n', A_H_cal);
end

fprintf('\nDone. Built both Tauchen and Rouwenhorst two-type baselines.\n');
