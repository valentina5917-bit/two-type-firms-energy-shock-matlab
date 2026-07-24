%German Industry in the Hopenhayn Model: the 2022 energy shock
%Two-firm-type crisis simulation using observed monthly price paths

clear; clc;

code_dir = fileparts(mfilename('fullpath'));
root = fileparts(code_dir);
addpath(code_dir, fullfile(code_dir, 'helpers'));
cd(root);

%Loading German steady state calibration
baseline_file = fullfile(root, 'data', 'baselines', ...
    'German_Baseline_TwoTypes_TAUCHEN.mat');
baseline = load(baseline_file);

required_fields = { ...
    'alpha','beta','x','xE','gamma_type','A_type','lambda_type', ...
    'rho_theta','sigma_theta', ...
    'W','Pe','L_supply','n_theta','theta_grid','Q_incumbent','Q_entrant', ...
    'V_type','exit_policy','mu','mu_operate','exit_rate_type','exit_rate_agg', ...
    'M_total','M_entrants_type','y_policy','n_policy','e_policy','pi_policy'};
missing_fields = setdiff(required_fields, fieldnames(baseline));
if ~isempty(missing_fields)
    error('Two-type baseline is missing: %s', strjoin(missing_fields, ', '));
end

alpha = baseline.alpha;
beta = baseline.beta;
x = baseline.x;
xE = baseline.xE;
gamma_type = baseline.gamma_type(:)';
A_type = baseline.A_type(:)';
lambda_type = baseline.lambda_type(:)';
W = baseline.W;
Pe = baseline.Pe;
L_supply = baseline.L_supply;
n_theta = baseline.n_theta;
theta_grid = baseline.theta_grid(:);
Q_incumbent = baseline.Q_incumbent;
Q_entrant = baseline.Q_entrant(:);
V_type = baseline.V_type;
exit_policy = baseline.exit_policy;
mu = baseline.mu;
mu_operate = baseline.mu_operate;
exit_rate_type = baseline.exit_rate_type(:);
exit_rate_agg = baseline.exit_rate_agg;
M_total = baseline.M_total;
M_entrants_type = baseline.M_entrants_type(:);
baseline_y_policy = baseline.y_policy;
baseline_n_policy = baseline.n_policy;
baseline_e_policy = baseline.e_policy;
baseline_pi_policy = baseline.pi_policy;

n_types = numel(gamma_type);
if n_types ~= 2
    error('Crisis_Simulation requires exactly two firm types.');
end
if ~isequal(size(mu), [n_theta n_types]) || ~isequal(size(V_type), [n_theta n_types])
    error('Two-type baseline arrays have incompatible dimensions.');
end

[~, idx_low] = min(gamma_type);
[~, idx_high] = max(gamma_type);

fprintf('\nLoaded two-type baseline: %s\n', baseline_file);
fprintf('Low-energy gamma  = %.8f\n', gamma_type(idx_low));
fprintf('High-energy gamma = %.8f\n', gamma_type(idx_high));
for s = 1:n_types
    eta_s = 1 - alpha - gamma_type(s);
    if eta_s <= 0
        error('Type %d has non-positive eta_s = 1-alpha-gamma_s.', s);
    end
    fprintf('Type %d eta = %.8f, output elasticity wrt Pe = %.8f\n', ...
        s, eta_s, -gamma_type(s)/eta_s);
end

%% Monthly calibratiom
%These conversions are local to the transition paths illustrative exercise. The annual baseline 
%file and the yearly IRFs remain unchanged.

periods_per_year = 12;
period_flow_scale = 1/periods_per_year;

beta_annual = beta;
x_annual = x;
L_supply_annual = L_supply;
rho_theta_annual = baseline.rho_theta;
sigma_theta_annual = baseline.sigma_theta;
Q_incumbent_annual = Q_incumbent;
Q_entrant_annual = Q_entrant;
V_type_annual = V_type;
exit_policy_annual = exit_policy;
mu_annual = mu;
mu_operate_annual = mu_operate;
exit_rate_type_annual = exit_rate_type;
exit_rate_agg_annual = exit_rate_agg;
M_entrants_type_annual = M_entrants_type;

beta = beta_annual^(1/periods_per_year);
x_month_flow = x_annual*period_flow_scale;
x = x_month_flow;
L_supply = L_supply_annual*period_flow_scale;

rho_theta = rho_theta_annual^(1/periods_per_year);
sigma_theta = sigma_theta_annual ...
    * sqrt((1-rho_theta^2)/(1-rho_theta_annual^2));
[z_month,Q_incumbent] = tauchen( ...
    n_theta,0,rho_theta,sigma_theta,4);
theta_grid = exp(z_month(:));
Q_entrant = stationary_markov_distribution(Q_incumbent);

monthly_grid_error = max(abs( ...
    log(theta_grid) - log(baseline.theta_grid(:))));
if monthly_grid_error > 1e-8
    warning('Monthly Tauchen grid differs from the annual grid by %.3e.', ...
        monthly_grid_error);
end

%Solving the baseline Bellman equations using monthly flows and monthly
%transitions. Entry remains a one-time cost, so xE is not divided by 12.

V_type = zeros(n_theta,n_types);
exit_policy = zeros(n_theta,n_types);
baseline_y_policy = zeros(n_theta,n_types);
baseline_n_policy = zeros(n_theta,n_types);
baseline_e_policy = zeros(n_theta,n_types);
baseline_pi_policy = zeros(n_theta,n_types);
monthly_baseline_iterations = zeros(n_types,1);
monthly_baseline_residual = zeros(n_types,1);

for s = 1:n_types
    gamma_s = gamma_type(s);
    A_s = A_type(s);
    eta_s = 1-alpha-gamma_s;

    y_base = period_flow_scale*(A_s*theta_grid).^(1/eta_s) ...
        .* (alpha/W)^(alpha/eta_s) ...
        .* (gamma_s/Pe)^(gamma_s/eta_s);
    n_base = alpha*y_base/W;
    e_base = gamma_s*y_base/Pe;
    pi_base = y_base-W*n_base-Pe*e_base;

    baseline_y_policy(:,s) = y_base;
    baseline_n_policy(:,s) = n_base;
    baseline_e_policy(:,s) = e_base;
    baseline_pi_policy(:,s) = pi_base;

end

target_type_mass = M_total*sum(mu_annual,1);
target_exit_rate_monthly = ...
    1-(1-exit_rate_agg_annual)^(1/periods_per_year);
monthly_fixed_cost_scale = calibrate_monthly_fixed_cost_scale( ...
    baseline_pi_policy,Q_incumbent,beta,W,x_month_flow, ...
    V_type_annual,Q_entrant,target_type_mass,target_exit_rate_monthly);
x = x_month_flow*monthly_fixed_cost_scale;

for s = 1:n_types
    [V_type(:,s),exit_policy(:,s),monthly_baseline_iterations(s), ...
        monthly_baseline_residual(s)] = solve_stationary_value( ...
        baseline_pi_policy(:,s),Q_incumbent,beta,W,x,V_type_annual(:,s));
end

% Reconstruct a monthly stationary distribution while preserving each
% type's calibrated mass. The implied monthly entry flow replaces monthly
% exits and therefore does not compound an annual entry flow twelve times.
monthly_pre_mass = zeros(n_theta,n_types);
monthly_operating_mass = zeros(n_theta,n_types);
M_entrants_type = zeros(n_types,1);

for s = 1:n_types
    operate_s = 1-exit_policy(:,s);
    survivor_transition = Q_incumbent'*diag(operate_s);
    unit_entry_mass = (eye(n_theta)-survivor_transition)\Q_entrant;
    if any(unit_entry_mass < -1e-10) || sum(unit_entry_mass) <= 0
        error('Could not construct a nonnegative monthly distribution for type %d.',s);
    end
    unit_entry_mass = max(unit_entry_mass,0);
    M_entrants_type(s) = target_type_mass(s)/sum(unit_entry_mass);
    monthly_pre_mass(:,s) = M_entrants_type(s)*unit_entry_mass;
    monthly_operating_mass(:,s) = monthly_pre_mass(:,s).*operate_s;
end

mu = monthly_pre_mass/M_total;
mu_operate = monthly_operating_mass/M_total;
exit_rate_type = sum(monthly_pre_mass.*exit_policy,1)' ...
    ./ sum(monthly_pre_mass,1)';
exit_rate_agg = sum(monthly_pre_mass(:).*exit_policy(:)) ...
    / sum(monthly_pre_mass(:));

exit_rate_type_annualized = 1-(1-exit_rate_type).^periods_per_year;
exit_rate_agg_annualized = 1-(1-exit_rate_agg)^periods_per_year;
monthly_stationarity_error = max(abs( ...
    Q_incumbent'*monthly_operating_mass ...
    + Q_entrant*M_entrants_type' - monthly_pre_mass),[],'all');

fprintf('\nMonthly crisis calibration:\n');
fprintf('beta annual / monthly          = %.8f / %.8f\n',beta_annual,beta);
fprintf('rho theta annual / monthly     = %.8f / %.8f\n', ...
    rho_theta_annual,rho_theta);
fprintf('monthly fixed-cost scale       = %.8f\n',monthly_fixed_cost_scale);
fprintf('monthly stationarity error     = %.3e\n',monthly_stationarity_error);
fprintf('aggregate exit monthly         = %.6f%%\n',100*exit_rate_agg);
fprintf('aggregate exit annualized      = %.6f%% (annual target %.6f%%)\n', ...
    100*exit_rate_agg_annualized,100*exit_rate_agg_annual);
fprintf('low exit annualized            = %.6f%% (annual target %.6f%%)\n', ...
    100*exit_rate_type_annualized(idx_low),100*exit_rate_type_annual(idx_low));
fprintf('high exit annualized           = %.6f%% (annual target %.6f%%)\n', ...
    100*exit_rate_type_annualized(idx_high),100*exit_rate_type_annual(idx_high));

%% Empirical data for energy-price and labor cost indices

observed_data_dir = fullfile(root, 'data', 'observed');
raw_pe_data = readmatrix(fullfile(observed_data_dir, 'Pet.csv'), ...
    'NumHeaderLines', 1);
raw_w_data = readmatrix(fullfile(observed_data_dir, 'wt.csv'), ...
    'NumHeaderLines', 1);

if size(raw_pe_data,2) < 2 || size(raw_w_data,2) < 2
    error('Pet.csv and wt.csv must each contain dates and one numeric series.');
end

Pe_full = raw_pe_data(:,2);
w_full = raw_w_data(:,2);

if numel(Pe_full) ~= numel(w_full)
    error('Energy-price and wage paths must have the same number of observations.');
end
if any(~isfinite(Pe_full)) || any(~isfinite(w_full)) ...
        || any(Pe_full <= 0) || any(w_full <= 0)
    error('Energy-price and wage paths must be finite and strictly positive.');
end

data_start_date = datetime(2017,1,1);
full_timeline = data_start_date + calmonths(0:numel(Pe_full)-1);
transition_start_date = datetime(2021,1,1);
transition_start_idx = find(full_timeline >= transition_start_date,1,'first');
if isempty(transition_start_idx)
    error('The transition start date lies outside the observed sample.');
end

%Retaining December 2020 only as the lag for January 2021 output growth

if transition_start_idx == 1
    error('A pre-transition observation is required to calculate initial growth.');
end
growth_lag_date = full_timeline(transition_start_idx-1);
Pe_growth_lag = Pe_full(transition_start_idx-1);
W_growth_lag = w_full(transition_start_idx-1);

Pe_observed = Pe_full(transition_start_idx:end);
w_observed = w_full(transition_start_idx:end);
Pe_series = Pe_observed;
w_series = w_observed;
timeline = full_timeline(transition_start_idx:end);
T = numel(Pe_series);

sc = 1;
% 1 = observed data
% 2 = observed path to the peak, then permanent peak energy price and rigid wage
% 3 = observed path to the peak, then mean reversion to the pre-crisis economy

export_compact_figure = true;
export_input_figure = true;

chi_common = 1.00;
chi_type = 0.00;
use_margin_exit_response = true;
exit_margin_scale = 0.75;
exit_margin_gain = 1.00;
exit_expectation_horizon_years = 50;
exit_expectation_horizon = periods_per_year*exit_expectation_horizon_years;
exit_expectation_rho_annual = 0.80;
exit_expectation_rho = exit_expectation_rho_annual^(1/periods_per_year);

%The observed crisis is treated as unanticipated before 2021. This keeps
% pre-crisis exits at their calibrated monthly steady-state rates.

crisis_start_date = transition_start_date;
crisis_start_idx = find(timeline >= crisis_start_date,1,'first');
if isempty(crisis_start_idx)
    error('The crisis start date lies outside the observed sample.');
end

vcheck_base_policy = zeros(n_theta,n_types);
for s = 1:n_types
    vcheck_base_policy(:,s) = -W*x + beta*(Q_incumbent*V_type(:,s));
end

switch sc
    case 1
        scenario_name = 'Actual data from 2021';
        scenario_suffix = 'actual_data_from_2021';

    case 2
        scenario_name = 'Permanent crisis and rigid wages';
        scenario_suffix = 'permanent_crisis';
        [~, peak_idx] = max(Pe_series);
        Pe_series(peak_idx:end) = Pe_series(peak_idx);
        w_series(:) = W;

    case 3
        scenario_name = 'Temporary crisis with recovery';
        scenario_suffix = 'temporary_recovery';
        [~, peak_idx] = max(Pe_series);
        rho_pe = 0.90^(1/periods_per_year);
        rho_w = 0.95^(1/periods_per_year);
        for t = peak_idx+1:T
            Pe_series(t) = Pe + rho_pe*(Pe_series(t-1) - Pe);
            w_series(t) = W + rho_w*(w_series(t-1) - W);
        end

    otherwise
        error('Invalid scenario. Select sc = 1, 2, or 3.');
end

disp(scenario_name);

Pe_terminal = Pe_series(end);
W_terminal = w_series(end);

%% Type specific terminal values

%The observed-data terminal environment is assumed to persist after the
%last observation. Each type therefore receives its own terminal value.

V_terminal = zeros(n_theta,n_types);
terminal_exit_policy = zeros(n_theta,n_types);
terminal_iterations = zeros(n_types,1);

policy_max_iter = 100;
bellman_tol = 1e-10;
terminal_start = tic;

for s = 1:n_types
    gamma_s = gamma_type(s);
    A_s = A_type(s);
    eta_s = 1 - alpha - gamma_s;
    shock_gap_terminal = max(Pe_terminal/Pe - 1,0);
    gamma_rel = gamma_s/mean(gamma_type);
    x_eff_terminal = x*(1 + chi_common*shock_gap_terminal ...
        + chi_type*shock_gap_terminal*gamma_rel);

    y_terminal = period_flow_scale*(A_s*theta_grid).^(1/eta_s) ...
        .* (alpha/W_terminal)^(alpha/eta_s) ...
        .* (gamma_s/Pe_terminal)^(gamma_s/eta_s);
    n_terminal = alpha*y_terminal/W_terminal;
    e_terminal = gamma_s*y_terminal/Pe_terminal;
    pi_terminal = y_terminal - W_terminal*n_terminal - Pe_terminal*e_terminal;

    V_guess = V_type(:,s);
    continue_policy = ...
        (-W_terminal*x_eff_terminal + beta*(Q_incumbent*V_guess)) > 0;
    policy_converged = false;

    for it = 1:policy_max_iter
        continuation_matrix = diag(double(continue_policy));
        policy_system = eye(n_theta) ...
            - beta*continuation_matrix*Q_incumbent;
        policy_rhs = pi_terminal ...
            - W_terminal*x_eff_terminal*double(continue_policy);
        V_candidate = policy_system\policy_rhs;

        new_continue_policy = ...
            (-W_terminal*x_eff_terminal + beta*(Q_incumbent*V_candidate)) > 0;
        if isequal(new_continue_policy,continue_policy)
            policy_converged = true;
            break;
        end
        continue_policy = new_continue_policy;
    end

    if ~policy_converged
        error('Terminal policy iteration did not converge for type %d.', s);
    end

    terminal_continuation = ...
        -W_terminal*x_eff_terminal + beta*(Q_incumbent*V_candidate);
    bellman_residual = max(abs(V_candidate ...
        - (pi_terminal + max(terminal_continuation,0))));
    if bellman_residual > bellman_tol
        error('Terminal Bellman residual is %.3e for type %d.', ...
            bellman_residual,s);
    end

    V_terminal(:,s) = V_candidate;
    terminal_iterations(s) = it;
    terminal_exit_policy(:,s) = double(terminal_continuation <= 0);
end

terminal_time = toc(terminal_start);
fprintf('Terminal prices: Pe = %.6f, W = %.6f\n', Pe_terminal, W_terminal);
fprintf('Terminal policy iterations: low = %d, high = %d\n', ...
    terminal_iterations(idx_low), terminal_iterations(idx_high));
fprintf('Terminal solution time: %.3f seconds\n',terminal_time);

%% Backward induction

disp('Starting backward induction...');
tic;

V_path = zeros(n_theta,T,n_types);
exit_policy_path = zeros(n_theta,T,n_types);
exit_prob_path = zeros(n_theta,T,n_types);
vcheck_path = zeros(n_theta,T,n_types);
y_policy_path = zeros(n_theta,T,n_types);
n_policy_path = zeros(n_theta,T,n_types);
e_policy_path = zeros(n_theta,T,n_types);
pi_policy_path = zeros(n_theta,T,n_types);

V_tomorrow = V_terminal;

for t = T:-1:1
    Pe_today = Pe_series(t);
    W_today = w_series(t);

    for s = 1:n_types
        gamma_s = gamma_type(s);
        A_s = A_type(s);
        eta_s = 1 - alpha - gamma_s;

        y_today = period_flow_scale*(A_s*theta_grid).^(1/eta_s) ...
            .* (alpha/W_today)^(alpha/eta_s) ...
            .* (gamma_s/Pe_today)^(gamma_s/eta_s);
        n_today = alpha*y_today/W_today;
        e_today = gamma_s*y_today/Pe_today;
        pi_today = y_today - W_today*n_today - Pe_today*e_today;

        y_policy_path(:,t,s) = y_today;
        n_policy_path(:,t,s) = n_today;
        e_policy_path(:,t,s) = e_today;
        pi_policy_path(:,t,s) = pi_today;

        shock_gap = max(Pe_today/Pe - 1,0);
        gamma_rel = gamma_s/mean(gamma_type);
        x_eff = x*(1 + chi_common*shock_gap ...
            + chi_type*shock_gap*gamma_rel);

        continuation = ...
            -W_today*x_eff + beta*(Q_incumbent*V_tomorrow(:,s));
        vcheck_path(:,t,s) = continuation;
        exit_policy_path(:,t,s) = double(continuation <= 0);
        V_path(:,t,s) = pi_today + max(continuation,0);
    end

    V_tomorrow = reshape(V_path(:,t,:), n_theta, n_types);
end

backward_time = toc;
fprintf('Backward induction complete in %.3f seconds.\n', backward_time);

%Firms expect that month's price disturbance to decay at the monthly
%equivalent of annual persistence rho = 0.80.

exit_mapping = struct( ...
    'alpha',alpha,'beta',beta,'x',x,'W',W,'Pe',Pe, ...
    'theta_grid',theta_grid,'Q_incumbent',Q_incumbent, ...
    'V_type',V_type,'exit_policy',exit_policy, ...
    'gamma_type',gamma_type,'A_type',A_type, ...
    'chi_common',chi_common,'chi_type',chi_type, ...
    'use_margin_exit_response',use_margin_exit_response, ...
    'exit_margin_scale',exit_margin_scale, ...
    'exit_margin_gain',exit_margin_gain, ...
    'vcheck_base_policy',vcheck_base_policy, ...
    'horizon',exit_expectation_horizon,'rho',exit_expectation_rho, ...
    'flow_scale',period_flow_scale);

reference_pre_mass = M_total*mu;
reference_type_mass = sum(reference_pre_mass,1);
reference_total_mass = sum(reference_type_mass);
exit_vcheck_impact_path = zeros(n_theta,T,n_types);
exit_rate_transmission_type_path = zeros(T,n_types);
exit_rate_transmission_agg_path = zeros(T,1);

for s = 1:n_types
    [exit_prob_s,vcheck_impact_s] = ...
        irf_consistent_exit_probabilities(Pe_series,s,exit_mapping);
    exit_prob_path(:,:,s) = exit_prob_s;
    exit_vcheck_impact_path(:,:,s) = vcheck_impact_s;

    if crisis_start_idx > 1
        pre_crisis_idx = 1:crisis_start_idx-1;
        exit_prob_path(:,pre_crisis_idx,s) = ...
            repmat(exit_policy(:,s),1,numel(pre_crisis_idx));
        exit_vcheck_impact_path(:,pre_crisis_idx,s) = ...
            repmat(vcheck_base_policy(:,s),1,numel(pre_crisis_idx));
    end

    exit_rate_transmission_type_path(:,s) = ...
        (reference_pre_mass(:,s)'*exit_prob_path(:,:,s) ...
        / reference_type_mass(s))';
    exit_rate_transmission_agg_path = exit_rate_transmission_agg_path ...
        + (reference_pre_mass(:,s)'*exit_prob_path(:,:,s))' ...
        / reference_total_mass;
end

%% Forward simulation

disp('Starting forward simulation...');
tic;

m_current = M_total*mu;
entry_flow_type = M_entrants_type;

mass_pre_total_path = zeros(T,1);
mass_pre_type_path = zeros(T,n_types);
mass_active_total_path = zeros(T,1);
mass_active_type_path = zeros(T,n_types);
mass_next_total_path = zeros(T,1);

exit_flow_total_path = zeros(T,1);
exit_flow_type_path = zeros(T,n_types);
exit_rate_agg_path = zeros(T,1);
exit_rate_type_path = zeros(T,n_types);

Y_agg_path = zeros(T,1);
N_agg_path = zeros(T,1);
E_agg_path = zeros(T,1);
Pi_agg_path = zeros(T,1);
Y_type_path = zeros(T,n_types);
N_type_path = zeros(T,n_types);
E_type_path = zeros(T,n_types);
Pi_type_path = zeros(T,n_types);

FC_agg_path = zeros(T,1);
C_path = zeros(T,1);
goods_gap = zeros(T,1);

for t = 1:T
    Pe_today = Pe_series(t);
    W_today = w_series(t);

    mass_pre_total_path(t) = sum(m_current(:));
    mass_pre_type_path(t,:) = sum(m_current,1);

    m_operate = zeros(n_theta,n_types);
    m_next = zeros(n_theta,n_types);

    for s = 1:n_types
        exit_prob_t = exit_prob_path(:,t,s);
        operate = 1 - exit_prob_t;
        m_operate(:,s) = m_current(:,s).*operate;

        mass_active_type_path(t,s) = sum(m_operate(:,s));
        exit_flow_type_path(t,s) = sum(m_current(:,s).*exit_prob_t);

        if mass_pre_type_path(t,s) > 0
            exit_rate_type_path(t,s) = exit_flow_type_path(t,s)/mass_pre_type_path(t,s);
        end

        Y_type_path(t,s) = sum(y_policy_path(:,t,s).*m_operate(:,s));
        N_type_path(t,s) = sum(n_policy_path(:,t,s).*m_operate(:,s));
        E_type_path(t,s) = sum(e_policy_path(:,t,s).*m_operate(:,s));
        Pi_type_path(t,s) = sum(pi_policy_path(:,t,s).*m_operate(:,s));

        m_next(:,s) = Q_incumbent'*m_operate(:,s) ...
            + entry_flow_type(s)*Q_entrant;
    end

    mass_active_total_path(t) = sum(mass_active_type_path(t,:));
    exit_flow_total_path(t) = sum(exit_flow_type_path(t,:));
    if mass_pre_total_path(t) > 0
        exit_rate_agg_path(t) = exit_flow_total_path(t)/mass_pre_total_path(t);
    end

    Y_agg_path(t) = sum(Y_type_path(t,:));
    N_agg_path(t) = sum(N_type_path(t,:));
    E_agg_path(t) = sum(E_type_path(t,:));
    Pi_agg_path(t) = sum(Pi_type_path(t,:));

    FC_agg_path(t) = W_today*(mass_active_total_path(t)*x + sum(entry_flow_type)*xE);
    C_path(t) = Pe_today*E_agg_path(t) + W_today*L_supply ...
        + Pi_agg_path(t) - FC_agg_path(t);
    goods_gap(t) = C_path(t) - Y_agg_path(t);

    m_next(m_next < 0) = 0;
    mass_next_total_path(t) = sum(m_next(:));
    m_current = m_next;
end

forward_time = toc;
fprintf('Forward simulation complete in %.3f seconds.\n', forward_time);

%% Baseline normalization
baseline_operating_mass = M_total*mu_operate;
baseline_active_type = sum(baseline_operating_mass,1);
baseline_active_total = sum(baseline_active_type);

Y0_type = sum(baseline_y_policy.*baseline_operating_mass,1);
N0_type = sum(baseline_n_policy.*baseline_operating_mass,1);
E0_type = sum(baseline_e_policy.*baseline_operating_mass,1);
Pi0_type = sum(baseline_pi_policy.*baseline_operating_mass,1);

Y0_agg = sum(Y0_type);
N0_agg = sum(N0_type);
E0_agg = sum(E0_type);
Pi0_agg = sum(Pi0_type);

Pe_index = 100*Pe_series/Pe;
W_index = 100*w_series/W;
exit_rate_agg_pct = 100*exit_rate_agg_path;
exit_rate_type_pct = 100*exit_rate_type_path;
exit_rate_transmission_agg_pct = 100*exit_rate_transmission_agg_path;
exit_rate_transmission_type_pct = 100*exit_rate_transmission_type_path;

%Base-100 normalization for comparison with the observed price path

if exit_rate_agg_path(1) <= 0 || any(exit_rate_type_path(1,:) <= 0)
    error('Initial exit rates must be positive for base-100 normalization.');
end
exit_rate_agg_index = 100*exit_rate_agg_path/exit_rate_agg_path(1);
exit_rate_type_index = zeros(T,n_types);
for s = 1:n_types
    exit_rate_type_index(:,s) = ...
        100*exit_rate_type_path(:,s)/exit_rate_type_path(1,s);
end

mass_active_total_index = 100*mass_active_total_path/baseline_active_total;
mass_active_type_index = zeros(T,n_types);
Y_type_index = zeros(T,n_types);
N_type_index = zeros(T,n_types);
E_type_index = zeros(T,n_types);

for s = 1:n_types
    mass_active_type_index(:,s) = 100*mass_active_type_path(:,s)/baseline_active_type(s);
    Y_type_index(:,s) = 100*Y_type_path(:,s)/Y0_type(s);
    N_type_index(:,s) = 100*N_type_path(:,s)/N0_type(s);
    E_type_index(:,s) = 100*E_type_path(:,s)/E0_type(s);
end

Y_index = 100*Y_agg_path/Y0_agg;
Y_dev_pct = Y_index - 100;
Y_type_dev_pct = Y_type_index - 100;
N_index = 100*N_agg_path/N0_agg;
E_index = 100*E_agg_path/E0_agg;

all_core_series = [exit_rate_agg_pct, exit_rate_type_pct, ...
    exit_rate_agg_index, exit_rate_type_index, ...
    mass_active_total_index, mass_active_type_index, ...
    Y_index, Y_type_index, E_index, E_type_index];
if any(~isfinite(all_core_series(:)))
    error('The crisis simulation generated non-finite output.');
end

%% Diagnostics

[peak_Pe_index, peak_Pe_idx] = max(Pe_index);
[peak_exit_agg, peak_exit_idx] = max(exit_rate_agg_pct);
exit_diff_pct = exit_rate_type_pct(:,idx_high) - exit_rate_type_pct(:,idx_low);
[peak_exit_diff, peak_exit_diff_idx] = max(exit_diff_pct);
[peak_exit_transmission, peak_exit_transmission_idx] = ...
    max(exit_rate_transmission_agg_pct);
exit_transmission_diff_pct = ...
    exit_rate_transmission_type_pct(:,idx_high) ...
    - exit_rate_transmission_type_pct(:,idx_low);
[peak_exit_transmission_diff, peak_exit_transmission_diff_idx] = ...
    max(exit_transmission_diff_pct);

fprintf('\n====================================================\n');
fprintf('TWO-TYPE CRISIS SIMULATION: %s\n', upper(scenario_name));
fprintf('====================================================\n');
fprintf('Observations                    = %d months\n', T);
fprintf('Peak energy-price index         = %.4f in %s\n', ...
    peak_Pe_index, datestr(timeline(peak_Pe_idx), 'mmm-yyyy'));
fprintf('Baseline monthly aggregate exit = %.6f%%\n', 100*exit_rate_agg);
fprintf('Exit mechanism                   = anchored marginal response\n');
fprintf('Crisis information date          = %s\n', ...
    datestr(crisis_start_date,'mmm-yyyy'));
fprintf('chi common / type               = %.3f / %.3f\n', chi_common,chi_type);
fprintf('Peak fixed-composition monthly exit = %.6f%% in %s\n', ...
    peak_exit_transmission, ...
    datestr(timeline(peak_exit_transmission_idx), 'mmm-yyyy'));
fprintf('Peak realized monthly exit      = %.6f%% in %s\n', ...
    peak_exit_agg, datestr(timeline(peak_exit_idx), 'mmm-yyyy'));
fprintf('Peak fixed-composition H-L gap  = %.6f p.p. in %s\n', ...
    peak_exit_transmission_diff, ...
    datestr(timeline(peak_exit_transmission_diff_idx), 'mmm-yyyy'));
fprintf('Exit index at peak Pe: aggregate=%.3f, low=%.3f, high=%.3f\n', ...
    exit_rate_agg_index(peak_Pe_idx), exit_rate_type_index(peak_Pe_idx,idx_low), ...
    exit_rate_type_index(peak_Pe_idx,idx_high));
fprintf('Fixed-composition monthly exit at peak Pe: agg=%.3f%%, low=%.3f%%, high=%.3f%%\n', ...
    exit_rate_transmission_agg_pct(peak_Pe_idx), ...
    exit_rate_transmission_type_pct(peak_Pe_idx,idx_low), ...
    exit_rate_transmission_type_pct(peak_Pe_idx,idx_high));
fprintf('Realized monthly exit at peak Pe: agg=%.3f%%, low=%.3f%%, high=%.3f%%\n', ...
    exit_rate_agg_pct(peak_Pe_idx),exit_rate_type_pct(peak_Pe_idx,idx_low), ...
    exit_rate_type_pct(peak_Pe_idx,idx_high));
fprintf('Initial aggregate output        = %.6f\n', Y_agg_path(1));
fprintf('Initial low output              = %.6f\n', Y_type_path(1,idx_low));
fprintf('Initial high output             = %.6f\n', Y_type_path(1,idx_high));
fprintf('Minimum aggregate output        = %.6f\n', min(Y_agg_path));
fprintf('Minimum low output              = %.6f\n', min(Y_type_path(:,idx_low)));
fprintf('Minimum high output             = %.6f\n', min(Y_type_path(:,idx_high)));
fprintf('Maximum absolute goods gap      = %.6e\n', max(abs(goods_gap)));
fprintf('====================================================\n');

%% PLOTS

window_size = min(6,T);

Y_growth_lag_type = zeros(1,n_types);
for s = 1:n_types
    gamma_s = gamma_type(s);
    eta_s = 1 - alpha - gamma_s;
    y_growth_lag = period_flow_scale*(A_type(s)*theta_grid).^(1/eta_s) ...
        .* (alpha/W_growth_lag)^(alpha/eta_s) ...
        .* (gamma_s/Pe_growth_lag)^(gamma_s/eta_s);
    Y_growth_lag_type(s) = ...
        sum(y_growth_lag.*baseline_operating_mass(:,s));
end
Y_growth_lag_agg = sum(Y_growth_lag_type);

Y_agg_growth_pct = zeros(T,1);
Y_type_growth_pct = zeros(T,n_types);
Y_agg_growth_pct(1) = 100*(Y_agg_path(1)/Y_growth_lag_agg - 1);
Y_type_growth_pct(1,:) = ...
    100*(Y_type_path(1,:)./Y_growth_lag_type - 1);
Y_agg_growth_pct(2:end) = 100*(Y_agg_path(2:end)./Y_agg_path(1:end-1) - 1);
Y_type_growth_pct(2:end,:) = ...
    100*(Y_type_path(2:end,:)./Y_type_path(1:end-1,:) - 1);
smoothed_exit_type = smoothdata(exit_rate_transmission_type_pct,1,'movmean',window_size);
smoothed_Y_type = smoothdata(Y_type_growth_pct,1,'movmean',window_size);

start_date = timeline(1);
custom_ticks = start_date:calmonths(6):timeline(end);

color_low = [0.00 0.20 0.70];
color_high = [0.70 0.05 0.10];
base_line = [0.45 0.45 0.45];

fig_crisis = figure('Color','w', ...
    'Name','Two-Type Crisis Paths', ...
    'Position',[160 40 980 820]);

tl_crisis = tiledlayout(fig_crisis,2,1, ...
    'TileSpacing','compact','Padding','compact');
tl_crisis.Units = 'normalized';
tl_crisis.Position = [0.09 0.22 0.88 0.75];

%Six-month moving averages of fixed-composition monthly exit paths

ax_exit = nexttile(tl_crisis);
h_exit_low = plot(ax_exit,timeline,smoothed_exit_type(:,idx_low), ...
    'Color',color_low,'LineWidth',1.2);
hold(ax_exit,'on');
h_exit_high = plot(ax_exit,timeline,smoothed_exit_type(:,idx_high), ...
    'Color',color_high,'LineWidth',1.2);
ylabel(ax_exit,'Monthly Exit Rate (%)');
all_exit_rates = exit_rate_transmission_type_pct(:);
ylim(ax_exit,[0, 1.10*max(all_exit_rates) + 0.1]);

%Six-month moving averages of month-to-month output growth by firm type
ax_output = nexttile(tl_crisis);
h_Y_low = plot(ax_output,timeline,smoothed_Y_type(:,idx_low), ...
    'Color',color_low,'LineWidth',1.2);
hold(ax_output,'on');
h_Y_high = plot(ax_output,timeline,smoothed_Y_type(:,idx_high), ...
    'Color',color_high,'LineWidth',1.2);
yline(ax_output,0,'--','Color',base_line,'LineWidth',1.2);
ylabel(ax_output,'Monthly Output Growth (%)');
ylim(ax_output,[-2 1.5]);
xlabel(ax_output,'Date');

annotation(fig_crisis,'rectangle',[0.15 0.025 0.70 0.060], ...
    'Color',[0 0 0],'LineWidth',1.0);
annotation(fig_crisis,'line',[0.18 0.22],[0.055 0.055], ...
    'Color',color_high,'LineWidth',1.2);
annotation(fig_crisis,'textbox',[0.23 0.035 0.275 0.040], ...
    'String','High energy-intensive firms','EdgeColor','none', ...
    'FontName','Times New Roman','FontSize',11, ...
    'VerticalAlignment','middle');
annotation(fig_crisis,'line',[0.52 0.56],[0.055 0.055], ...
    'Color',color_low,'LineWidth',1.2);
annotation(fig_crisis,'textbox',[0.57 0.035 0.265 0.040], ...
    'String','Low energy-intensive firms','EdgeColor','none', ...
    'FontName','Times New Roman','FontSize',11, ...
    'VerticalAlignment','middle');

crisis_axes = [ax_exit ax_output];
for aa = 1:numel(crisis_axes)
    grid(crisis_axes(aa),'off');
    box(crisis_axes(aa),'off');
    crisis_axes(aa).XAxisLocation = 'bottom';
    crisis_axes(aa).YAxisLocation = 'left';
    crisis_axes(aa).FontName = 'Times New Roman';
    crisis_axes(aa).FontSize = 13;
    crisis_axes(aa).LineWidth = 1.1;
    crisis_axes(aa).TickDir = 'out';
    crisis_axes(aa).XLim = [timeline(1) timeline(end)];
    crisis_axes(aa).XTick = custom_ticks;
    crisis_axes(aa).XAxis.TickLabelFormat = 'MM/yy';
    crisis_axes(aa).XAxis.TickLabelRotation = 90;
    try
        crisis_axes(aa).Toolbar.Visible = 'off';
    catch
    end
end
ax_exit.YLabel.FontSize = 13;
ax_output.YLabel.FontSize = 13;
ax_output.XLabel.FontSize = 13;

function p_exit = marginal_exit_probability( ...
        vcheck,vcheck_base,hard_exit,scale,gain)
    p_base = smooth_exit_probability(vcheck_base,scale);
    p_shock = smooth_exit_probability(vcheck,scale);
    extra_exit = gain*max(p_shock-p_base,0);
    p_exit = hard_exit+(1-hard_exit).*min(extra_exit,1);
end

function p_exit = smooth_exit_probability(vcheck,scale)
    z = max(min(vcheck./scale,40),-40);
    p_exit = 1./(1+exp(z));
end

function [p_exit,vcheck_impact] = ...
        irf_consistent_exit_probabilities(Pe_impacts,s,p)
    gamma_s = p.gamma_type(s);
    A_s = p.A_type(s);
    eta_s = 1-p.alpha-gamma_s;
    gamma_rel = gamma_s/mean(p.gamma_type);

    Pe_impacts = reshape(Pe_impacts,1,[]);
    n_paths = numel(Pe_impacts);
    log_impact = log(Pe_impacts/p.Pe);
    V_next = repmat(p.V_type(:,s),1,n_paths);
    vcheck_impact = repmat(p.vcheck_base_policy(:,s),1,n_paths);
    productivity_term = (A_s*p.theta_grid).^(1/eta_s) ...
        *(p.alpha/p.W)^(p.alpha/eta_s);

    for h = p.horizon:-1:1
        Pe_h = p.Pe*exp((p.rho^(h-1))*log_impact);
        energy_price_term = (gamma_s./Pe_h).^(gamma_s/eta_s);
        y_h = p.flow_scale*(productivity_term*energy_price_term);
        n_h = p.alpha*y_h/p.W;
        e_h = gamma_s*y_h./Pe_h;
        pi_h = y_h-p.W*n_h-e_h.*Pe_h;

        shock_gap = max(Pe_h/p.Pe-1,0);
        x_eff = p.x*(1+p.chi_common*shock_gap ...
            +p.chi_type*shock_gap*gamma_rel);
        vcheck_h = -p.W*x_eff+p.beta*(p.Q_incumbent*V_next);
        V_current = pi_h+max(vcheck_h,0);

        if h == 1
            vcheck_impact = vcheck_h;
        end
        V_next = V_current;
    end

    if p.use_margin_exit_response
        p_exit = marginal_exit_probability( ...
            vcheck_impact,p.vcheck_base_policy(:,s),p.exit_policy(:,s), ...
            p.exit_margin_scale,p.exit_margin_gain);
    else
        p_exit = double(vcheck_impact <= 0);
    end
end

function q = stationary_markov_distribution(Q)
    [vectors,values] = eig(Q');
    [~,idx] = min(abs(diag(values)-1));
    q = real(vectors(:,idx));
    if sum(q) < 0
        q = -q;
    end
    q = max(q,0);
    q = q/sum(q);
    residual = norm(Q'*q-q,inf);
    if residual > 1e-10
        error('Monthly entrant distribution residual is %.3e.',residual);
    end
end

function scale = calibrate_monthly_fixed_cost_scale( ...
        pi_policy,Q,beta_period,W,x_month_flow,V_guess,Q_entrant, ...
        target_type_mass,target_exit_rate)
    lower = 0.25;
    upper = 4.00;
    lower_rate = monthly_exit_rate_for_cost( ...
        lower,pi_policy,Q,beta_period,W,x_month_flow,V_guess, ...
        Q_entrant,target_type_mass);
    upper_rate = monthly_exit_rate_for_cost( ...
        upper,pi_policy,Q,beta_period,W,x_month_flow,V_guess, ...
        Q_entrant,target_type_mass);

    while upper_rate < target_exit_rate && upper < 64
        upper = 2*upper;
        upper_rate = monthly_exit_rate_for_cost( ...
            upper,pi_policy,Q,beta_period,W,x_month_flow,V_guess, ...
            Q_entrant,target_type_mass);
    end
    if lower_rate > target_exit_rate || upper_rate < target_exit_rate
        error('Could not bracket the monthly fixed-cost calibration.');
    end

    candidate_scales = zeros(38,1);
    candidate_rates = zeros(38,1);
    candidate_scales(1:2) = [lower; upper];
    candidate_rates(1:2) = [lower_rate; upper_rate];

    for it = 1:36
        middle = 0.5*(lower+upper);
        middle_rate = monthly_exit_rate_for_cost( ...
            middle,pi_policy,Q,beta_period,W,x_month_flow,V_guess, ...
            Q_entrant,target_type_mass);
        candidate_scales(it+2) = middle;
        candidate_rates(it+2) = middle_rate;

        if middle_rate < target_exit_rate
            lower = middle;
        else
            upper = middle;
        end
    end

    [~,best_idx] = min(abs(candidate_rates-target_exit_rate));
    scale = candidate_scales(best_idx);
end

function aggregate_rate = monthly_exit_rate_for_cost( ...
        scale,pi_policy,Q,beta_period,W,x_month_flow,V_guess, ...
        Q_entrant,target_type_mass)
    n = size(Q,1);
    n_types_local = size(pi_policy,2);
    type_rate = zeros(n_types_local,1);

    for s = 1:n_types_local
        [~,policy_s] = solve_stationary_value( ...
            pi_policy(:,s),Q,beta_period,W,x_month_flow*scale,V_guess(:,s));
        if ~any(policy_s)
            type_rate(s) = 0;
            continue;
        end

        survivor_transition = Q'*diag(1-policy_s);
        unit_mass = (eye(n)-survivor_transition)\Q_entrant;
        unit_mass = max(real(unit_mass),0);
        type_rate(s) = sum(unit_mass.*policy_s)/sum(unit_mass);
    end

    aggregate_rate = sum(target_type_mass(:).*type_rate) ...
        / sum(target_type_mass);
end

function [V,exit_policy,it,residual] = solve_stationary_value( ...
        pi_flow,Q,beta_period,W,x_period,V_guess)
    n = numel(pi_flow);
    continue_policy = ...
        (-W*x_period+beta_period*(Q*V_guess)) > 0;
    converged = false;

    for it = 1:100
        continuation_matrix = diag(double(continue_policy));
        system_matrix = eye(n)-beta_period*continuation_matrix*Q;
        rhs = pi_flow-W*x_period*double(continue_policy);
        V_candidate = system_matrix\rhs;
        new_continue_policy = ...
            (-W*x_period+beta_period*(Q*V_candidate)) > 0;
        if isequal(new_continue_policy,continue_policy)
            converged = true;
            break;
        end
        continue_policy = new_continue_policy;
    end

    if ~converged
        error('Monthly baseline policy iteration did not converge.');
    end

    continuation = -W*x_period+beta_period*(Q*V_candidate);
    V = pi_flow+max(continuation,0);
    residual = max(abs(V-V_candidate));
    exit_policy = double(continuation <= 0);
end
