%German Industry in the Hopenhayn Model
%Two-type one-shock simulation with AR(1) recovery in energy prices

clear; clc;

code_dir = fileparts(mfilename('fullpath'));
root = fileparts(code_dir);
addpath(code_dir, fullfile(code_dir, 'helpers'));
cd(root);

%% Loading two-type steady state

baseline_file = fullfile(root, 'data', 'baselines', ...
    'German_Baseline_TwoTypes_TAUCHEN.mat');
load(baseline_file);

% baseline objects:
% V_type, mu, exit_policy, alpha, beta, x, xE,
% theta_grid, n_theta, Q_incumbent, Q_entrant,
% W, Pe, L_supply,
% gamma_L, gamma_H, gamma_type,
% lambda_H, lambda_type,
% A_L, A_H, A_type,
% exit_rate_type, exit_rate_agg,
% M_total, M_type, M_entrants_total, M_entrants_type

n_types = length(gamma_type);
[~, idx_high] = max(gamma_type);
idx_low = 3 - idx_high;

fprintf('\nLoaded two-type baseline.\n');
fprintf('gamma low  = %.8f\n', gamma_type(idx_low));
fprintf('gamma high = %.8f\n', gamma_type(idx_high));
fprintf('alpha      = %.8f\n', alpha);
fprintf('A_L        = %.8f\n', A_type(idx_low));
fprintf('A_H        = %.8f\n', A_type(idx_high));
fprintf('lambda low = %.8f\n', lambda_type(idx_low));
fprintf('lambda high= %.8f\n', lambda_type(idx_high));

for s = 1:n_types
    eta_s = 1 - alpha - gamma_type(s);
    fprintf('Type %d output elasticity wrt Pe = %.8f\n', s, -gamma_type(s)/eta_s);
end

%% Energy price shock

T = 50;                 % years
T_plot = 50;            

shock_size = 0.20;       % 20% initial increase in energy price

%Shock sizes for robustness
shock_size_robust_values = [0.10 0.20 0.40];

%Yearly persistence the shock
rho_pe = 0.80;

%Persistence values for robustness
rho_robust_values = [0.50 0.80 0.95];

%Composition firm for robustness
%NaN is replaced below by the exact calibrated baseline share
high_share_robust_values = [0.20 NaN 0.50];

%Fixed-cost wedge in survival condition
chi_common = 1.00;
chi_type   = 0.00;

% Extra marginal exit response for baseline survivors whose continuation
% values fall during the shock
use_margin_exit_response = true;
exit_margin_scale = 0.75;
exit_margin_gain  = 1.00;

%Reduced-form productivity/demand effect set to 0 for a pure energy-price shock
zeta_A = 0.0;

W_path = ones(T,1) * W;

logPe_path = zeros(T,1);
logPe_path(1) = log(1 + shock_size);

for t = 2:T
    logPe_path(t) = rho_pe * logPe_path(t-1);
end

Pe_path = Pe * exp(logPe_path);

baseline_exit_prob = exit_policy;
vcheck_base_policy = zeros(n_theta, n_types);

for s = 1:n_types
    vcheck_base_policy(:,s) = -W*x + beta * (Q_incumbent * V_type(:,s));
end

V_terminal = V_type;

%% Backward induction

disp('Starting backward induction...');
tic;

V_path           = zeros(n_theta, T, n_types);
exit_policy_path = zeros(n_theta, T, n_types); % 1 = exit, 0 = continue
exit_prob_path   = zeros(n_theta, T, n_types);
vcheck_path      = zeros(n_theta, T, n_types);

y_policy_path  = zeros(n_theta, T, n_types);
n_policy_path  = zeros(n_theta, T, n_types);
e_policy_path  = zeros(n_theta, T, n_types);
pi_policy_path = zeros(n_theta, T, n_types);


V_tomorrow = V_terminal;

for t = T:-1:1

    W_today  = W_path(t);
    Pe_today = Pe_path(t);

    for s = 1:n_types

        gamma_s = gamma_type(s);
        A_s     = A_type(s);

        shock_gap = max(Pe_today/Pe - 1, 0);

        A_today = A_s * (1 - zeta_A * shock_gap);
        A_today = max(A_today, 1e-6);

        eta_s = 1 - alpha - gamma_s;

        for i = 1:n_theta

            th = theta_grid(i);

            y_today = (A_today * th)^(1/eta_s) ...
                    * (alpha / W_today)^(alpha / eta_s) ...
                    * (gamma_s / Pe_today)^(gamma_s / eta_s);

            n_today = alpha * y_today / W_today;
            e_today = gamma_s * y_today / Pe_today;

            pi_today = y_today - W_today*n_today - Pe_today*e_today;

            y_policy_path(i,t,s)  = y_today;
            n_policy_path(i,t,s)  = n_today;
            e_policy_path(i,t,s)  = e_today;
            pi_policy_path(i,t,s) = pi_today;
        end

        EV = beta * (Q_incumbent * V_tomorrow(:,s));

        %Energy-linked fixed-cost wedge in continuation condition
        %High-energy firms receive a larger survival-cost shock
       
        shock_gap = max(Pe_today/Pe - 1, 0);
        gamma_rel = gamma_s / mean(gamma_type);

        x_eff = x * (1 + chi_common*shock_gap + chi_type*shock_gap*gamma_rel);

        %continuation decision
        vcheck = -W_today*x_eff + EV;
        vcheck_path(:,t,s) = vcheck;

        %exit rule for the value function, matching ss
        exit_policy_path(:,t,s) = double(vcheck <= 0);

        if use_margin_exit_response
            exit_prob_path(:,t,s) = marginal_exit_probability( ...
                vcheck, vcheck_base_policy(:,s), exit_policy(:,s), ...
                exit_margin_scale, exit_margin_gain);
        else
            exit_prob_path(:,t,s) = exit_policy_path(:,t,s);
        end

        vhat = max(vcheck, 0);

        %Value function
        V_path(:,t,s) = pi_policy_path(:,t,s) + vhat;
    end

    V_tomorrow = squeeze(V_path(:,t,:));
end

toc;
disp('Backward induction complete.');

fprintf('\nHard-exit diagnostic, impact Year t=1:\n');

for s = 1:n_types
    fprintf('Type %d: exit states = %d out of %d\n', ...
        s, sum(exit_policy_path(:,1,s) == 1), n_theta);
end

if use_margin_exit_response
    fprintf('\nAnchored marginal-exit accounting, impact Year t=1:\n');

    for s = 1:n_types
        margin_impact_exit = sum((M_total * mu(:,s)) .* exit_prob_path(:,1,s)) ...
            / sum(M_total * mu(:,s));
        fprintf('Type %d: expected exit rate = %.6f%%\n', ...
            s, 100*margin_impact_exit);
    end
end

%% Year 0 baseline inputs
m0_pre = M_total * mu;   %pre-exit stock at Year 0

%Baseline static policies at Pe = 1, W = 1
y0_policy  = zeros(n_theta, n_types);
n0_policy  = zeros(n_theta, n_types);
e0_policy  = zeros(n_theta, n_types);
pi0_policy = zeros(n_theta, n_types);

for s = 1:n_types

    gamma_s = gamma_type(s);
    A_s     = A_type(s);
    eta_s   = 1 - alpha - gamma_s;

    for i = 1:n_theta

        th = theta_grid(i);

        y0_policy(i,s) = (A_s * th)^(1/eta_s) ...
                       * (alpha / W)^(alpha / eta_s) ...
                       * (gamma_s / Pe)^(gamma_s / eta_s);

        n0_policy(i,s) = alpha * y0_policy(i,s) / W;
        e0_policy(i,s) = gamma_s * y0_policy(i,s) / Pe;

        pi0_policy(i,s) = y0_policy(i,s) ...
                        - W*n0_policy(i,s) ...
                        - Pe*e0_policy(i,s);
    end
end

%Year-0 operating firms from baseline exit rule
m0_operate = zeros(n_theta, n_types);

for s = 1:n_types
    m0_operate(:,s) = m0_pre(:,s) .* (1 - baseline_exit_prob(:,s));
end

M0_type_active = zeros(1,n_types);
Y0_type  = zeros(1,n_types);
N0_type  = zeros(1,n_types);
E0_type  = zeros(1,n_types);
Pi0_type = zeros(1,n_types);

for s = 1:n_types
    M0_type_active(s) = sum(m0_operate(:,s));
    Y0_type(s)  = sum(y0_policy(:,s)  .* m0_operate(:,s));
    N0_type(s)  = sum(n0_policy(:,s)  .* m0_operate(:,s));
    E0_type(s)  = sum(e0_policy(:,s)  .* m0_operate(:,s));
    Pi0_type(s) = sum(pi0_policy(:,s) .* m0_operate(:,s));
end

M0_active = sum(M0_type_active);

Y0_agg  = sum(Y0_type);
N0_agg  = sum(N0_type);
E0_agg  = sum(E0_type);
Pi0_agg = sum(Pi0_type);

%% Survival cut-off diagnostics
fprintf('\nSURVIVAL-CUTOFF DIAGNOSTICS AT IMPACT\n');

for s = 1:n_types

    gamma_s = gamma_type(s);
    gamma_rel = gamma_s / mean(gamma_type);

    Pe_today = Pe_path(1);
    W_today  = W_path(1);

    shock_gap = max(Pe_today/Pe - 1, 0);

    x_eff_impact = x * (1 + chi_common*shock_gap + chi_type*shock_gap*gamma_rel);

    EV_impact = beta * (Q_incumbent * squeeze(V_path(:,2,s)));
    vcheck_shock = -W_today*x_eff_impact + EV_impact;

    EV_base = beta * (Q_incumbent * V_type(:,s));
    vcheck_base = -W*x + EV_base;

    fprintf('\nType %d\n', s);
    fprintf('gamma_s                             = %.8f\n', gamma_s);
    fprintf('gamma_rel                           = %.8f\n', gamma_rel);
    fprintf('x_eff impact / x - 1                = %.8f%%\n', 100*(x_eff_impact/x - 1));
    fprintf('baseline exit states                = %d of %d\n', sum(vcheck_base <= 0), n_theta);
    fprintf('impact exit states                  = %d of %d\n', sum(vcheck_shock <= 0), n_theta);
    fprintf('newly exiting states                = %d of %d\n', sum((vcheck_base > 0) & (vcheck_shock <= 0)), n_theta);

    if use_margin_exit_response
        exit_prob_base = exit_policy(:,s);
        exit_prob_shock = marginal_exit_probability( ...
            vcheck_shock, vcheck_base_policy(:,s), exit_policy(:,s), ...
            exit_margin_scale, exit_margin_gain);

        margin_base_rate = sum(m0_pre(:,s) .* exit_prob_base) / sum(m0_pre(:,s));
        margin_shock_rate = sum(m0_pre(:,s) .* exit_prob_shock) / sum(m0_pre(:,s));

        fprintf('margin baseline exit                = %.6f%%\n', 100*margin_base_rate);
        fprintf('margin impact exit                  = %.6f%%\n', 100*margin_shock_rate);
    end

    if any(vcheck_base > 0)
        fprintf('min baseline vcheck among survivors = %.8e\n', min(vcheck_base(vcheck_base > 0)));
    else
        fprintf('min baseline vcheck among survivors = NaN\n');
    end

    if any(vcheck_shock > 0)
        fprintf('min shock vcheck among survivors    = %.8e\n', min(vcheck_shock(vcheck_shock > 0)));
    else
        fprintf('min shock vcheck among survivors    = NaN\n');
    end
end

%% Forward simulation

disp('Starting forward simulation...');
tic;

m_current = m0_pre;

entry_flow_type_ss  = M_entrants_type(:);
entry_flow_total_ss = sum(entry_flow_type_ss);

mass_pre_total_path    = zeros(T+1,1);
mass_pre_type_path     = zeros(T+1,n_types);

mass_active_total_path = zeros(T+1,1);
mass_active_type_path  = zeros(T+1,n_types);

exit_flow_total_path = zeros(T,1);
exit_flow_type_path  = zeros(T,n_types);

entry_flow_total_path = zeros(T,1);
entry_flow_type_path  = zeros(T,n_types);

exit_rate_agg_path  = zeros(T+1,1);
exit_rate_type_path = zeros(T+1,n_types);

Y_agg_path  = zeros(T+1,1);
N_agg_path  = zeros(T+1,1);
E_agg_path  = zeros(T+1,1);
Pi_agg_path = zeros(T+1,1);
FC_agg_path = zeros(T+1,1);

Y_type_path  = zeros(T+1,n_types);
N_type_path  = zeros(T+1,n_types);
E_type_path  = zeros(T+1,n_types);
Pi_type_path = zeros(T+1,n_types);

%Year 0
mass_pre_total_path(1) = sum(m0_pre(:));

for s = 1:n_types
    mass_pre_type_path(1,s) = sum(m0_pre(:,s));
end

mass_active_total_path(1)  = M0_active;
mass_active_type_path(1,:) = M0_type_active;

exit_rate_agg_path(1)    = exit_rate_agg;
exit_rate_type_path(1,:) = exit_rate_type(:)';

Y_agg_path(1)  = Y0_agg;
N_agg_path(1)  = N0_agg;
E_agg_path(1)  = E0_agg;
Pi_agg_path(1) = Pi0_agg;

%Baseline fixed costs: operating firms + baseline entrants
FC_agg_path(1) = W * (M0_active*x + entry_flow_total_ss*xE);

Y_type_path(1,:)  = Y0_type;
N_type_path(1,:)  = N0_type;
E_type_path(1,:)  = E0_type;
Pi_type_path(1,:) = Pi0_type;

for t = 1:T

    %Pre-exit stock at start of ear
    mass_pre_total_path(t+1) = sum(m_current(:));

    for s = 1:n_types
        mass_pre_type_path(t+1,s) = sum(m_current(:,s));
    end

    m_operate = zeros(n_theta, n_types);
    m_next    = zeros(n_theta, n_types);

    for s = 1:n_types

        exit_prob_t = exit_prob_path(:,t,s);
        operate = 1 - exit_prob_t;
        m_operate(:,s) = m_current(:,s) .* operate;

        %Only operating firms produce
        Y_type_path(t+1,s)  = sum(y_policy_path(:,t,s)  .* m_operate(:,s));
        N_type_path(t+1,s)  = sum(n_policy_path(:,t,s)  .* m_operate(:,s));
        E_type_path(t+1,s)  = sum(e_policy_path(:,t,s)  .* m_operate(:,s));
        Pi_type_path(t+1,s) = sum(pi_policy_path(:,t,s) .* m_operate(:,s));

        %Exit flows and rates relative to pre-exit stock
        exit_flow_type_path(t,s) = sum(m_current(:,s) .* exit_prob_t);

        if mass_pre_type_path(t+1,s) > 0
            exit_rate_type_path(t+1,s) = exit_flow_type_path(t,s) / mass_pre_type_path(t+1,s);
        else
            exit_rate_type_path(t+1,s) = 0;
        end

        %Entry closure: replace exiting firms of the same type
        entry_flow_type_t = exit_flow_type_path(t,s);
        entry_flow_type_path(t,s) = entry_flow_type_t;

        %continuing incumbents + replacement entrants
        m_next(:,s) = Q_incumbent' * m_operate(:,s) + entry_flow_type_t * Q_entrant;
    end

    entry_flow_total_path(t) = sum(entry_flow_type_path(t,:));

    mass_active_total_path(t+1) = sum(m_operate(:));

    for s = 1:n_types
        mass_active_type_path(t+1,s) = sum(m_operate(:,s));
    end

    Y_agg_path(t+1)  = sum(Y_type_path(t+1,:));
    N_agg_path(t+1)  = sum(N_type_path(t+1,:));
    E_agg_path(t+1)  = sum(E_type_path(t+1,:));
    Pi_agg_path(t+1) = sum(Pi_type_path(t+1,:));

    exit_flow_total_path(t) = sum(exit_flow_type_path(t,:));

    if mass_pre_total_path(t+1) > 0
        exit_rate_agg_path(t+1) = exit_flow_total_path(t) / mass_pre_total_path(t+1);
    else
        exit_rate_agg_path(t+1) = 0;
    end

    %operating firms + time-varying entrants
    operating_mass_total = sum(m_operate(:));
    FC_agg_path(t+1) = W_path(t) * (operating_mass_total*x + entry_flow_total_path(t)*xE);

    m_next(m_next < 0) = 0;

    m_current = m_next;
end

toc;
disp('Forward simulation complete.');

%% Indexes and IRFs

years = 0:T;
plot_end = min(T_plot, T);
plot_idx = 1:(plot_end+1);

Pe_path_full = [Pe; Pe_path];
Pe_index = 100 * Pe_path_full / Pe;

exit_rate_agg_pct  = 100 * exit_rate_agg_path;
exit_rate_type_pct = 100 * exit_rate_type_path;

mass_active_total_index = 100 * mass_active_total_path / M0_active;

mass_active_type_index = zeros(T+1,n_types);

for s = 1:n_types
    mass_active_type_index(:,s) = 100 * mass_active_type_path(:,s) / M0_type_active(s);
end

Y_index = 100 * Y_agg_path / Y0_agg;
N_index = 100 * N_agg_path / N0_agg;
E_index = 100 * E_agg_path / E0_agg;

%Type-specific indexes
Y_type_index = zeros(T+1,n_types);
N_type_index = zeros(T+1,n_types);
E_type_index = zeros(T+1,n_types);

for s = 1:n_types
    Y_type_index(:,s) = 100 * Y_type_path(:,s) / Y0_type(s);
    N_type_index(:,s) = 100 * N_type_path(:,s) / N0_type(s);
    E_type_index(:,s) = 100 * E_type_path(:,s) / E0_type(s);
end

%% Smoothing
smooth_window = 1;

%Raw series
exit_rate_agg_pct_plot       = exit_rate_agg_pct;
exit_rate_type_pct_plot      = exit_rate_type_pct;

mass_active_total_index_plot = mass_active_total_index;
mass_active_type_index_plot  = mass_active_type_index;

Y_index_plot      = Y_index;
Y_type_index_plot = Y_type_index;

E_index_plot      = E_index;
E_type_index_plot = E_type_index;

exit_rate_agg_pct_plot(2:end) = smoothdata(exit_rate_agg_pct(2:end), ...
    'gaussian', smooth_window);

mass_active_total_index_plot(2:end) = smoothdata(mass_active_total_index(2:end), ...
    'gaussian', smooth_window);

Y_index_plot(2:end) = smoothdata(Y_index(2:end), ...
    'gaussian', smooth_window);

E_index_plot(2:end) = smoothdata(E_index(2:end), ...
    'gaussian', smooth_window);

for s = 1:n_types

    exit_rate_type_pct_plot(2:end,s) = smoothdata(exit_rate_type_pct(2:end,s), ...
        'gaussian', smooth_window);

    mass_active_type_index_plot(2:end,s) = smoothdata(mass_active_type_index(2:end,s), ...
        'gaussian', smooth_window);

    Y_type_index_plot(2:end,s) = smoothdata(Y_type_index(2:end,s), ...
        'gaussian', smooth_window);

    E_type_index_plot(2:end,s) = smoothdata(E_type_index(2:end,s), ...
        'gaussian', smooth_window);
end

%% Summary

if use_margin_exit_response
    exit_rule_label = 'ANCHORED MARGINAL EXIT';
    exit_rule_suffix = 'marginexit';
else
    exit_rule_label = 'HARD EXIT';
    exit_rule_suffix = 'hardexit';
end

fprintf('\n====================================================\n');
fprintf('TWO-TYPE ONE-SHOCK SIMULATION - %s, YEARLY MODEL\n', exit_rule_label);
fprintf('====================================================\n');
fprintf('shock size in Pe              = %+6.2f%%\n', 100*shock_size);
fprintf('rho_pe yearly                 = %.3f\n', rho_pe);
fprintf('chi_common fixed-cost wedge   = %.3f\n', chi_common);
fprintf('chi_type fixed-cost wedge     = %.3f\n', chi_type);
fprintf('zeta_A productivity wedge     = %.3f\n', zeta_A);
fprintf('margin exit response          = %d\n', use_margin_exit_response);
fprintf('exit margin scale             = %.3f\n', exit_margin_scale);
fprintf('exit margin gain              = %.3f\n', exit_margin_gain);
fprintf('entry rule                    = time-varying entrants by type\n');
fprintf('mass anchor normalization     = OFF\n');
fprintf('--------------------------------------------\n');
fprintf('baseline aggregate exit       = %.6f%%\n', 100*exit_rate_agg);
fprintf('impact-year agg exit          = %.6f%%\n', exit_rate_agg_pct(2));
fprintf('peak agg exit                 = %.6f%%\n', max(exit_rate_agg_pct));
fprintf('--------------------------------------------\n');
fprintf('baseline low-type exit        = %.6f%%\n', 100*exit_rate_type(idx_low));
fprintf('baseline high-type exit       = %.6f%%\n', 100*exit_rate_type(idx_high));
fprintf('impact-year low-type exit     = %.6f%%\n', exit_rate_type_pct(2,idx_low));
fprintf('impact-year high-type exit    = %.6f%%\n', exit_rate_type_pct(2,idx_high));
fprintf('--------------------------------------------\n');
fprintf('active mass index after 1y    = %.6f\n', mass_active_total_index(min(2,T+1)));
fprintf('output index after 1y         = %.6f\n', Y_index(min(2,T+1)));
fprintf('energy index after 1y         = %.6f\n', E_index(min(2,T+1)));
fprintf('active mass index after 5y    = %.6f\n', mass_active_total_index(min(6,T+1)));
fprintf('output index after 5y         = %.6f\n', Y_index(min(6,T+1)));
fprintf('energy index after 5y         = %.6f\n', E_index(min(6,T+1)));
fprintf('impact-year aggregate Y       = %.6f\n', Y_index(2));
fprintf('impact-year low-type Y        = %.6f\n', Y_type_index(2,idx_low));
fprintf('impact-year high-type Y       = %.6f\n', Y_type_index(2,idx_high));
fprintf('====================================================\n');

%% Checks
disp('Year-0 index checks, should be 100:');
disp([mass_active_total_index(1), Y_index(1), N_index(1), E_index(1)]);

disp('Year-0 mass active by type, should be 100:');
disp(mass_active_type_index(1,:));

disp('Year-0 output by type, should be 100:');
disp(Y_type_index(1,:));

disp('Year-0 energy by type, should be 100:');
disp(E_type_index(1,:));

fprintf('\nRaw impact output by type:\n');
fprintf('Low-energy Y_index Year 1  = %.8f\n', Y_type_index(2,idx_low));
fprintf('High-energy Y_index Year 1 = %.8f\n', Y_type_index(2,idx_high));
fprintf('Aggregate Y_index Year 1   = %.8f\n', Y_index(2));

fprintf('\nRaw impact energy use by type:\n');
fprintf('Low-energy E_index Year 1  = %.8f\n', E_type_index(2,idx_low));
fprintf('High-energy E_index Year 1 = %.8f\n', E_type_index(2,idx_high));
fprintf('Aggregate E_index Year 1   = %.8f\n', E_index(2));

%% PLOTS

plot_end = min(T_plot, T);

shock_plot_idx = 2:(plot_end+1);     
tgrid_raw      = 0:(plot_end-1);     
tgrid_fine     = linspace(0, plot_end-1, 8*plot_end);

% Interpolate for cleaner plotting only
Pe_index_fine = interp1(tgrid_raw, Pe_index(shock_plot_idx), ...
    tgrid_fine, 'pchip');

Exit_agg_plot_fine = interp1(tgrid_raw, exit_rate_agg_pct_plot(shock_plot_idx), ...
    tgrid_fine, 'pchip');

Exit_low_plot_fine = interp1(tgrid_raw, exit_rate_type_pct_plot(shock_plot_idx,idx_low), ...
    tgrid_fine, 'pchip');

Exit_high_plot_fine = interp1(tgrid_raw, exit_rate_type_pct_plot(shock_plot_idx,idx_high), ...
    tgrid_fine, 'pchip');

Mass_total_plot_fine = interp1(tgrid_raw, mass_active_total_index_plot(shock_plot_idx), ...
    tgrid_fine, 'pchip');

Mass_low_plot_fine = interp1(tgrid_raw, mass_active_type_index_plot(shock_plot_idx,idx_low), ...
    tgrid_fine, 'pchip');

Mass_high_plot_fine = interp1(tgrid_raw, mass_active_type_index_plot(shock_plot_idx,idx_high), ...
    tgrid_fine, 'pchip');

Y_index_plot_fine = interp1(tgrid_raw, Y_index_plot(shock_plot_idx), ...
    tgrid_fine, 'pchip');

Y_low_plot_fine = interp1(tgrid_raw, Y_type_index_plot(shock_plot_idx,idx_low), ...
    tgrid_fine, 'pchip');

Y_high_plot_fine = interp1(tgrid_raw, Y_type_index_plot(shock_plot_idx,idx_high), ...
    tgrid_fine, 'pchip');

E_index_plot_fine = interp1(tgrid_raw, E_index_plot(shock_plot_idx), ...
    tgrid_fine, 'pchip');

E_low_plot_fine = interp1(tgrid_raw, E_type_index_plot(shock_plot_idx,idx_low), ...
    tgrid_fine, 'pchip');

E_high_plot_fine = interp1(tgrid_raw, E_type_index_plot(shock_plot_idx,idx_high), ...
    tgrid_fine, 'pchip');


Pe_dev_fine = Pe_index_fine - 100;

Exit_agg_dev_fine = Exit_agg_plot_fine - 100*exit_rate_agg;
Exit_low_dev_fine = Exit_low_plot_fine - 100*exit_rate_type(idx_low);
Exit_high_dev_fine = Exit_high_plot_fine - 100*exit_rate_type(idx_high);

Mass_total_dev_fine = Mass_total_plot_fine - 100;
Mass_low_dev_fine = Mass_low_plot_fine - 100;
Mass_high_dev_fine = Mass_high_plot_fine - 100;

Y_dev_fine = Y_index_plot_fine - 100;
Y_low_dev_fine = Y_low_plot_fine - 100;
Y_high_dev_fine = Y_high_plot_fine - 100;

E_dev_fine = E_index_plot_fine - 100;
E_low_dev_fine = E_low_plot_fine - 100;
E_high_dev_fine = E_high_plot_fine - 100;

fig_symbolic = figure('Color','w', ...
       'Name',sprintf('IRFs to a %.0f%% energy-price shock',100*shock_size), ...
       'Position',[160 60 980 780]);

tl_symbolic = tiledlayout(fig_symbolic,3,3,'TileSpacing','compact','Padding','compact');
tl_symbolic.Units = 'normalized';
tl_symbolic.Position = [0.065 0.215 0.91 0.735];

line_agg = {'k-', 'LineWidth', 1.8};
color_low = [0.00 0.20 0.70];
color_high = [0.70 0.05 0.10];
line_low = {'-', 'Color', color_low, 'LineWidth', 1.8};
line_high = {'-', 'Color', color_high, 'LineWidth', 1.8};
base_line = [0.45 0.45 0.45];

nexttile;
plot(tgrid_fine, Pe_dev_fine, line_agg{:}); hold on;
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$P_e$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Exit_agg_dev_fine, line_agg{:}); hold on;
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$\delta$','Interpreter','latex');
ylabel('p.p. from ss');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Exit_low_dev_fine, line_low{:}); hold on;
plot(tgrid_fine, Exit_high_dev_fine, line_high{:});
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$\delta_s$','Interpreter','latex');
ylabel('p.p. from ss');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Mass_total_dev_fine, line_agg{:}); hold on;
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$M$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Mass_low_dev_fine, line_low{:}); hold on;
plot(tgrid_fine, Mass_high_dev_fine, line_high{:});
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$M_s$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Y_dev_fine, line_agg{:}); hold on;
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$Y$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, Y_low_dev_fine, line_low{:}); hold on;
plot(tgrid_fine, Y_high_dev_fine, line_high{:});
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$Y_s$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
xlabel('Years');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, E_dev_fine, line_agg{:}); hold on;
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$E$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
xlabel('Years');
format_symbolic_axis(plot_end);

nexttile;
plot(tgrid_fine, E_low_dev_fine, line_low{:}); hold on;
plot(tgrid_fine, E_high_dev_fine, line_high{:});
yline(0,'--','Color',base_line,'LineWidth',1.2);
title('$E_s$','Interpreter','latex');
ylabel('\% dev. from ss','Interpreter','latex');
xlabel('Years');
format_symbolic_axis(plot_end);

annotation(fig_symbolic, 'rectangle', [0.130 0.062 0.740 0.055], ...
    'Color', [0 0 0], 'LineWidth', 1.0);
annotation(fig_symbolic, 'line', [0.160 0.205], [0.090 0.090], ...
    'Color', 'k', 'LineWidth', 1.8);
annotation(fig_symbolic, 'textbox', [0.215 0.072 0.130 0.034], ...
    'String', 'Aggregate', ...
    'EdgeColor', 'none', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'VerticalAlignment', 'middle');
annotation(fig_symbolic, 'line', [0.385 0.430], [0.090 0.090], ...
    'Color', color_low, 'LineWidth', 1.8);
annotation(fig_symbolic, 'textbox', [0.440 0.072 0.175 0.034], ...
    'String', 'Low-energy firms', ...
    'EdgeColor', 'none', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'VerticalAlignment', 'middle');
annotation(fig_symbolic, 'line', [0.635 0.680], [0.090 0.090], ...
    'Color', color_high, 'LineWidth', 1.8);
annotation(fig_symbolic, 'textbox', [0.690 0.072 0.180 0.034], ...
    'String', 'High-energy firms', ...
    'EdgeColor', 'none', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'VerticalAlignment', 'middle');

caption_symbolic = sprintf('', ...
    100*shock_size);
annotation(fig_symbolic, 'textbox', [0.20 0.015 0.60 0.025], ...
    'String', caption_symbolic, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 13, ...
    'FontWeight', 'bold');

%% Robustness: shock persistence

rho_robust_values = unique(rho_robust_values(:)', 'stable');
n_rho_robust = numel(rho_robust_values);

rho_params = struct();
rho_params.T = T;
rho_params.n_theta = n_theta;
rho_params.n_types = n_types;
rho_params.shock_size = shock_size;
rho_params.rho_pe = rho_pe;
rho_params.W = W;
rho_params.Pe = Pe;
rho_params.x = x;
rho_params.beta = beta;
rho_params.alpha = alpha;
rho_params.theta_grid = theta_grid;
rho_params.gamma_type = gamma_type;
rho_params.A_type = A_type;
rho_params.V_type = V_type;
rho_params.exit_policy = exit_policy;
rho_params.exit_rate_type = exit_rate_type;
rho_params.M_total = M_total;
rho_params.mu = mu;
rho_params.Q_incumbent = Q_incumbent;
rho_params.Q_entrant = Q_entrant;
rho_params.idx_low = idx_low;
rho_params.idx_high = idx_high;
rho_params.chi_common = chi_common;
rho_params.chi_type = chi_type;
rho_params.zeta_A = zeta_A;
rho_params.use_margin_exit_response = use_margin_exit_response;
rho_params.exit_margin_scale = exit_margin_scale;
rho_params.exit_margin_gain = exit_margin_gain;

rho_robust_irfs = simulate_rho_persistence_irf(rho_robust_values(1), rho_params);

for rr = 2:n_rho_robust
    rho_robust_irfs(rr) = simulate_rho_persistence_irf(rho_robust_values(rr), rho_params);
end

rho_robust_tables = cell(n_rho_robust,1);

for rr = 1:n_rho_robust
    n_obs = numel(rho_robust_irfs(rr).years);
    rho_robust_tables{rr} = table( ...
        repmat(rho_robust_irfs(rr).rho_pe, n_obs, 1), ...
        rho_robust_irfs(rr).years(:), ...
        rho_robust_irfs(rr).Y_index(:), ...
        rho_robust_irfs(rr).Y_low_index(:), ...
        rho_robust_irfs(rr).Y_high_index(:), ...
        rho_robust_irfs(rr).Y_diff_high_minus_low_indexpts(:), ...
        rho_robust_irfs(rr).E_index(:), ...
        rho_robust_irfs(rr).E_low_index(:), ...
        rho_robust_irfs(rr).E_high_index(:), ...
        rho_robust_irfs(rr).exit_rate_low_pct(:), ...
        rho_robust_irfs(rr).exit_rate_high_pct(:), ...
        rho_robust_irfs(rr).exit_diff_high_minus_low_pctpts(:), ...
        rho_robust_irfs(rr).mass_active_low_index(:), ...
        rho_robust_irfs(rr).mass_active_high_index(:), ...
        rho_robust_irfs(rr).mass_active_diff_high_minus_low_index(:), ...
        'VariableNames', {'rho_pe','year','Y_index','Y_low_index','Y_high_index', ...
        'Y_diff_high_minus_low_indexpts','E_index','E_low_index','E_high_index', ...
        'exit_low_pct','exit_high_pct','exit_diff_high_minus_low_pctpts', ...
        'mass_active_low_index','mass_active_high_index','mass_active_diff_high_minus_low_index'});
end

rho_robustness_out = vertcat(rho_robust_tables{:});
rho_robust_colors = rho_persistence_colors(n_rho_robust);

fig_rho_robust = figure('Color','w', ...
    'Name','Robustness to Shock Persistence: Firm-Type IRFs', ...
    'Position',[120 40 980 860]);

tl_rho_robust = tiledlayout(fig_rho_robust,3,2,'TileSpacing','loose','Padding','loose');
tl_rho_robust.Units = 'normalized';
tl_rho_robust.Position = [0.075 0.220 0.90 0.665];

add_firm_type_column_headers(fig_rho_robust);

plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_low_index', 100, '% dev. from ss', '$Y_L$', plot_end, false);
plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_high_index', 100, '% dev. from ss', '$Y_H$', plot_end, false);
plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_rate_low_pct', 100*exit_rate_type(idx_low), 'p.p. from ss', '$\delta_L$', plot_end, false);
plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_rate_high_pct', 100*exit_rate_type(idx_high), 'p.p. from ss', '$\delta_H$', plot_end, false);
plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'mass_active_low_index', 100, '% dev. from ss', '$M_L$', plot_end, true);
plot_type_robustness_panel(rho_robust_irfs, rho_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'mass_active_high_index', 100, '% dev. from ss', '$M_H$', plot_end, true);

add_rho_persistence_legend(fig_rho_robust, rho_robust_values, rho_robust_colors);


%% Robustness: shock size

shock_size_robust_values = unique(shock_size_robust_values(:)', 'stable');
n_shock_size_robust = numel(shock_size_robust_values);

shock_size_params = rho_params;
shock_size_params.rho_pe = rho_pe;
shock_size_params.shock_size = shock_size_robust_values(1);
shock_size_robust_irfs = simulate_rho_persistence_irf(rho_pe, shock_size_params);

for ss = 2:n_shock_size_robust
    shock_size_params = rho_params;
    shock_size_params.rho_pe = rho_pe;
    shock_size_params.shock_size = shock_size_robust_values(ss);
    shock_size_robust_irfs(ss) = simulate_rho_persistence_irf(rho_pe, shock_size_params);
end

shock_size_robust_tables = cell(n_shock_size_robust,1);

for ss = 1:n_shock_size_robust
    n_obs = numel(shock_size_robust_irfs(ss).years);
    shock_size_robust_tables{ss} = table( ...
        repmat(shock_size_robust_irfs(ss).shock_size, n_obs, 1), ...
        shock_size_robust_irfs(ss).years(:), ...
        shock_size_robust_irfs(ss).Y_index(:), ...
        shock_size_robust_irfs(ss).Y_low_index(:), ...
        shock_size_robust_irfs(ss).Y_high_index(:), ...
        shock_size_robust_irfs(ss).Y_diff_high_minus_low_indexpts(:), ...
        shock_size_robust_irfs(ss).E_index(:), ...
        shock_size_robust_irfs(ss).E_low_index(:), ...
        shock_size_robust_irfs(ss).E_high_index(:), ...
        shock_size_robust_irfs(ss).exit_rate_low_pct(:), ...
        shock_size_robust_irfs(ss).exit_rate_high_pct(:), ...
        shock_size_robust_irfs(ss).exit_diff_high_minus_low_pctpts(:), ...
        shock_size_robust_irfs(ss).mass_active_low_index(:), ...
        shock_size_robust_irfs(ss).mass_active_high_index(:), ...
        shock_size_robust_irfs(ss).mass_active_diff_high_minus_low_index(:), ...
        'VariableNames', {'shock_size','year','Y_index','Y_low_index','Y_high_index', ...
        'Y_diff_high_minus_low_indexpts','E_index','E_low_index','E_high_index', ...
        'exit_low_pct','exit_high_pct','exit_diff_high_minus_low_pctpts', ...
        'mass_active_low_index','mass_active_high_index','mass_active_diff_high_minus_low_index'});
end

shock_size_robustness_out = vertcat(shock_size_robust_tables{:});
shock_size_robust_colors = rho_persistence_colors(n_shock_size_robust);

fig_shock_size_robust = figure('Color','w', ...
    'Name','Robustness to Shock Size: Firm-Type IRFs', ...
    'Position',[140 60 980 860]);

tl_shock_size_robust = tiledlayout(fig_shock_size_robust,3,2,'TileSpacing','loose','Padding','loose');
tl_shock_size_robust.Units = 'normalized';
tl_shock_size_robust.Position = [0.075 0.220 0.90 0.665];

add_firm_type_column_headers(fig_shock_size_robust);

plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_low_index', 100, '% dev. from ss', '$Y_L$', plot_end, false);
plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_high_index', 100, '% dev. from ss', '$Y_H$', plot_end, false);
plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_rate_low_pct', 100*exit_rate_type(idx_low), 'p.p. from ss', '$\delta_L$', plot_end, false);
plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_rate_high_pct', 100*exit_rate_type(idx_high), 'p.p. from ss', '$\delta_H$', plot_end, false);
plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'mass_active_low_index', 100, '% dev. from ss', '$M_L$', plot_end, true);
plot_type_robustness_panel(shock_size_robust_irfs, shock_size_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'mass_active_high_index', 100, '% dev. from ss', '$M_H$', plot_end, true);

add_shock_size_legend(fig_shock_size_robust, shock_size_robust_values, shock_size_robust_colors);


%% Robustness: firm-type composition

baseline_type_shares = sum(mu,1) / sum(mu(:));
baseline_high_share = baseline_type_shares(idx_high);
high_share_robust_values(isnan(high_share_robust_values)) = baseline_high_share;
high_share_robust_values = unique(max(min(high_share_robust_values(:)', 0.95), 0.05), 'stable');
n_high_share_robust = numel(high_share_robust_values);

share_params = rho_params;
share_params.rho_pe = rho_pe;
share_params.shock_size = shock_size;
share_params.mu = reweight_type_shares(mu, idx_low, idx_high, high_share_robust_values(1));
share_robust_irfs = simulate_rho_persistence_irf(rho_pe, share_params);

for hs = 2:n_high_share_robust
    share_params = rho_params;
    share_params.rho_pe = rho_pe;
    share_params.shock_size = shock_size;
    share_params.mu = reweight_type_shares(mu, idx_low, idx_high, high_share_robust_values(hs));
    share_robust_irfs(hs) = simulate_rho_persistence_irf(rho_pe, share_params);
end

share_robust_tables = cell(n_high_share_robust,1);

for hs = 1:n_high_share_robust
    n_obs = numel(share_robust_irfs(hs).years);
    share_robust_tables{hs} = table( ...
        repmat(share_robust_irfs(hs).high_type_share, n_obs, 1), ...
        share_robust_irfs(hs).years(:), ...
        share_robust_irfs(hs).Y_index(:), ...
        share_robust_irfs(hs).Y_low_index(:), ...
        share_robust_irfs(hs).Y_high_index(:), ...
        share_robust_irfs(hs).Y_diff_high_minus_low_indexpts(:), ...
        share_robust_irfs(hs).E_index(:), ...
        share_robust_irfs(hs).E_low_index(:), ...
        share_robust_irfs(hs).E_high_index(:), ...
        share_robust_irfs(hs).exit_rate_agg_pct(:), ...
        share_robust_irfs(hs).exit_rate_low_pct(:), ...
        share_robust_irfs(hs).exit_rate_high_pct(:), ...
        share_robust_irfs(hs).exit_diff_high_minus_low_pctpts(:), ...
        share_robust_irfs(hs).mass_active_index(:), ...
        share_robust_irfs(hs).mass_active_low_index(:), ...
        share_robust_irfs(hs).mass_active_high_index(:), ...
        share_robust_irfs(hs).mass_active_diff_high_minus_low_index(:), ...
        'VariableNames', {'high_type_share','year','Y_index','Y_low_index','Y_high_index', ...
        'Y_diff_high_minus_low_indexpts','E_index','E_low_index','E_high_index', ...
        'exit_agg_pct','exit_low_pct','exit_high_pct','exit_diff_high_minus_low_pctpts', ...
        'mass_active_index','mass_active_low_index','mass_active_high_index','mass_active_diff_high_minus_low_index'});
end

share_robustness_out = vertcat(share_robust_tables{:});
share_robust_colors = rho_persistence_colors(n_high_share_robust);

fig_share_robust = figure('Color','w', ...
    'Name','Robustness to Firm-Type Composition: Firm-Type IRFs', ...
    'Position',[90 15 1200 1040]);

tl_share_robust = tiledlayout(fig_share_robust,4,n_high_share_robust, ...
    'TileSpacing','compact','Padding','compact');
tl_share_robust.Units = 'normalized';
tl_share_robust.Position = [0.055 0.175 0.925 0.690];

column_left = tl_share_robust.Position(1);
column_width = tl_share_robust.Position(3) / n_high_share_robust;
[~, baseline_share_idx] = min(abs(high_share_robust_values - baseline_high_share));

for hs = 1:n_high_share_robust
    if hs == baseline_share_idx
        column_label = sprintf('$s_H = %.0f\\%%\\;\\mathrm{(baseline)}$', ...
            100*high_share_robust_values(hs));
    else
        column_label = sprintf('$s_H = %.0f\\%%$', 100*high_share_robust_values(hs));
    end

    annotation(fig_share_robust, 'textbox', ...
        [column_left+(hs-1)*column_width 0.928 column_width 0.032], ...
        'String', column_label, ...
        'Interpreter', 'latex', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'FontWeight', 'bold');
end

for hs = 1:n_high_share_robust
    plot_composition_scenario_panel(share_robust_irfs(hs), ...
        'Y_index', 'Y_low_index', 'Y_high_index', [100 100 100], ...
        '% dev. from ss', '$Y$', shock_plot_idx, tgrid_raw, tgrid_fine, plot_end, false);
end

for hs = 1:n_high_share_robust
    exit_baselines = [share_robust_irfs(hs).exit_rate_agg_pct(1), ...
        share_robust_irfs(hs).exit_rate_low_pct(1), ...
        share_robust_irfs(hs).exit_rate_high_pct(1)];
    plot_composition_scenario_panel(share_robust_irfs(hs), ...
        'exit_rate_agg_pct', 'exit_rate_low_pct', 'exit_rate_high_pct', exit_baselines, ...
        'p.p. from ss', '$\delta$', shock_plot_idx, tgrid_raw, tgrid_fine, plot_end, false);
end

for hs = 1:n_high_share_robust
    plot_composition_scenario_panel(share_robust_irfs(hs), ...
        'mass_active_index', 'mass_active_low_index', 'mass_active_high_index', [100 100 100], ...
        '% dev. from ss', '$M$', shock_plot_idx, tgrid_raw, tgrid_fine, plot_end, false);
end

for hs = 1:n_high_share_robust
    plot_composition_scenario_panel(share_robust_irfs(hs), ...
        'E_index', 'E_low_index', 'E_high_index', [100 100 100], ...
        '% dev. from ss', '$E$', shock_plot_idx, tgrid_raw, tgrid_fine, plot_end, true);
end

add_firm_type_legend(fig_share_robust);


%% Firm-type composition common normalization
% Each numerator is the shock-induced change relative to that scenario's
% own steady state. All scenarios use the calibrated composition's steady
% state as the denominator. The output gap uses calibrated aggregate output
% as its scale. The exit panel reports the high-minus-low difference in
% composition-weighted shock contributions.

common_norm_denominators = struct();
common_norm_denominators.high_type_share = baseline_high_share;
common_norm_denominators.Y = share_robust_irfs(baseline_share_idx).Y_ss;
common_norm_denominators.exit = share_robust_irfs(baseline_share_idx).exit_rate_agg_pct(1);
common_norm_denominators.M = share_robust_irfs(baseline_share_idx).mass_active_ss;
common_norm_denominators.E = share_robust_irfs(baseline_share_idx).E_ss;

share_common_irfs = struct([]);
share_common_tables = cell(n_high_share_robust,1);

for hs = 1:n_high_share_robust
    source_irf = share_robust_irfs(hs);

    share_common_irfs(hs).high_type_share = source_irf.high_type_share;
    share_common_irfs(hs).years = source_irf.years;
    share_common_irfs(hs).Y_common = 100 * source_irf.Y_ss ...
        .* (source_irf.Y_index/100 - 1) / common_norm_denominators.Y;
    share_common_irfs(hs).exit_common = 100 ...
        * (source_irf.exit_rate_agg_pct - source_irf.exit_rate_agg_pct(1)) ...
        / common_norm_denominators.exit;
    share_common_irfs(hs).M_common = 100 * source_irf.mass_active_ss ...
        .* (source_irf.mass_active_index/100 - 1) / common_norm_denominators.M;
    share_common_irfs(hs).E_common = 100 * source_irf.E_ss ...
        .* (source_irf.E_index/100 - 1) / common_norm_denominators.E;

    Y_high_change = source_irf.Y_high_ss .* (source_irf.Y_high_index/100 - 1);
    Y_low_change = source_irf.Y_low_ss .* (source_irf.Y_low_index/100 - 1);
    share_common_irfs(hs).Y_diff_common = 100 ...
        * (Y_high_change - Y_low_change) / common_norm_denominators.Y;

    exit_high_change = source_irf.exit_rate_high_pct - source_irf.exit_rate_high_pct(1);
    exit_low_change = source_irf.exit_rate_low_pct - source_irf.exit_rate_low_pct(1);
    high_share = source_irf.high_type_share;
    share_common_irfs(hs).exit_weighted_diff_common = 100 ...
        * (high_share*exit_high_change - (1-high_share)*exit_low_change) ...
        / common_norm_denominators.exit;

    common_fields = {'Y_common','exit_common','M_common','E_common', ...
        'Y_diff_common','exit_weighted_diff_common'};
    for ff = 1:numel(common_fields)
        share_common_irfs(hs).(common_fields{ff})(1) = 0;
    end

    n_obs = numel(source_irf.years);
    share_common_tables{hs} = table( ...
        repmat(source_irf.high_type_share, n_obs, 1), ...
        source_irf.years(:), ...
        share_common_irfs(hs).Y_common(:), ...
        share_common_irfs(hs).exit_common(:), ...
        share_common_irfs(hs).M_common(:), ...
        share_common_irfs(hs).E_common(:), ...
        share_common_irfs(hs).Y_diff_common(:), ...
        share_common_irfs(hs).exit_weighted_diff_common(:), ...
        'VariableNames', {'high_type_share','year','Y_common_pct','exit_common_pct', ...
        'M_common_pct','E_common_pct','Y_diff_high_minus_low_common_pct', ...
        'exit_weighted_high_minus_low_common_pct'});
end

share_common_normalization_out = vertcat(share_common_tables{:});

peak_share_labels = cell(n_high_share_robust,1);
peak_Y_common = zeros(n_high_share_robust,1);
peak_exit_common = zeros(n_high_share_robust,1);
peak_M_common = zeros(n_high_share_robust,1);
peak_E_common = zeros(n_high_share_robust,1);

for hs = 1:n_high_share_robust
    share_pct = 100 * share_common_irfs(hs).high_type_share;
    if hs == baseline_share_idx
        peak_share_labels{hs} = sprintf('%.0f\\%% (baseline)', share_pct);
    else
        peak_share_labels{hs} = sprintf('%.0f\\%%', share_pct);
    end

    peak_Y_common(hs) = signed_peak_response(share_common_irfs(hs).Y_common(shock_plot_idx));
    peak_exit_common(hs) = signed_peak_response(share_common_irfs(hs).exit_common(shock_plot_idx));
    peak_M_common(hs) = signed_peak_response(share_common_irfs(hs).M_common(shock_plot_idx));
    peak_E_common(hs) = signed_peak_response(share_common_irfs(hs).E_common(shock_plot_idx));
end

share_common_peak_table = table( ...
    reshape([share_common_irfs.high_type_share], [], 1), ...
    peak_share_labels, peak_Y_common, peak_exit_common, peak_M_common, peak_E_common, ...
    'VariableNames', {'high_type_share','share_label','peak_Y_pct', ...
    'peak_delta_pct','peak_M_pct','peak_E_pct'});

fig_share_common = figure('Color','w', ...
    'Name','Firm-Type Composition: Common-Normalized IRFs', ...
    'Position',[170 45 980 860]);

tl_share_common = tiledlayout(fig_share_common,3,2,'TileSpacing','loose','Padding','loose');
tl_share_common.Units = 'normalized';
tl_share_common.Position = [0.075 0.220 0.90 0.720];

common_norm_label = '% of baseline ss';
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_common', common_norm_label, '$Y$', plot_end, false);
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_common', common_norm_label, '$\delta$', plot_end, false);
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'M_common', common_norm_label, '$M$', plot_end, false);
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'E_common', common_norm_label, '$E$', plot_end, false);
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'Y_diff_common', common_norm_label, ...
    '$\Delta Y_H^{\mathrm{tot}} - \Delta Y_L^{\mathrm{tot}}$', plot_end, true);
plot_diff_robustness_panel(share_common_irfs, share_robust_colors, shock_plot_idx, tgrid_raw, tgrid_fine, ...
    'exit_weighted_diff_common', common_norm_label, ...
    '$s_H\Delta\delta_H - (1-s_H)\Delta\delta_L$', plot_end, true);

add_high_share_legend(fig_share_common, high_share_robust_values, share_robust_colors);


%% =========================================================
% 17. EXIT-DIFFERENTIAL DIAGNOSTICS
% ==========================================================
exit_diff_pct = exit_rate_type_pct(:,idx_high) - exit_rate_type_pct(:,idx_low);

[peak_diff, peak_diff_idx] = max(exit_diff_pct);

impact_diff = exit_diff_pct(2); 
avg_diff_5y = mean(exit_diff_pct(2:min(6,T+1)));

fprintf('\nEXIT-DIFFERENTIAL DIAGNOSTICS, HIGH - LOW\n');
fprintf('Impact differential, year 1:       %.6f p.p.\n', impact_diff);
fprintf('Peak differential:                 %.6f p.p. at year %d\n', peak_diff, peak_diff_idx-1);
fprintf('Average differential, years 1-5:   %.6f p.p.\n', avg_diff_5y);

if impact_diff <= 0
    warning('High-energy exit is NOT above low-energy exit on impact year.');
end

if peak_diff <= 0
    warning('High-energy exit is never above low-energy exit over the IRF horizon.');
end

%% Saving IRFs

irf_out = table(years(:), Pe_index(:), exit_rate_agg_pct(:), ...
    exit_rate_type_pct(:,idx_low), exit_rate_type_pct(:,idx_high), exit_diff_pct(:), ...
    mass_pre_total_path(:), mass_active_total_index(:), ...
    Y_index(:), E_index(:), ...
    Y_type_index(:,idx_low), Y_type_index(:,idx_high), ...
    E_type_index(:,idx_low), E_type_index(:,idx_high), ...
    'VariableNames', {'year','Pe_index','exit_agg_pct','exit_low_pct','exit_high_pct', ...
    'exit_diff_high_minus_low_pctpts','mass_pre_total','mass_active_index','Y_index','E_index', ...
    'Y_low_index','Y_high_index','E_low_index','E_high_index'});

results_folder = fullfile(root, 'results', 'generated');

if exist(results_folder, 'dir') ~= 7
    mkdir(results_folder);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_yearly_%s.csv', exit_rule_suffix, timestamp));

mat_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_yearly_%s.mat', exit_rule_suffix, timestamp));

symbolic_png_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_symbolic_%s.png', exit_rule_suffix, timestamp));

symbolic_pdf_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_symbolic_%s.pdf', exit_rule_suffix, timestamp));

rho_robust_csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_rho_robustness_%s.csv', exit_rule_suffix, timestamp));

rho_robust_mat_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_rho_robustness_%s.mat', exit_rule_suffix, timestamp));

rho_robust_png_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_rho_robustness_%s.png', exit_rule_suffix, timestamp));

rho_robust_pdf_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_rho_robustness_%s.pdf', exit_rule_suffix, timestamp));

shock_size_robust_csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_shock_size_robustness_%s.csv', exit_rule_suffix, timestamp));

shock_size_robust_mat_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_shock_size_robustness_%s.mat', exit_rule_suffix, timestamp));

shock_size_robust_png_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_shock_size_robustness_%s.png', exit_rule_suffix, timestamp));

shock_size_robust_pdf_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_shock_size_robustness_%s.pdf', exit_rule_suffix, timestamp));

share_robust_csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_robustness_%s.csv', exit_rule_suffix, timestamp));

share_robust_mat_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_robustness_%s.mat', exit_rule_suffix, timestamp));

share_robust_png_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_robustness_%s.png', exit_rule_suffix, timestamp));

share_robust_pdf_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_robustness_%s.pdf', exit_rule_suffix, timestamp));

share_common_csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_common_normalization_%s.csv', exit_rule_suffix, timestamp));

share_common_mat_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_common_normalization_%s.mat', exit_rule_suffix, timestamp));

share_common_png_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_common_normalization_%s.png', exit_rule_suffix, timestamp));

share_common_pdf_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_common_normalization_%s.pdf', exit_rule_suffix, timestamp));

share_common_peak_csv_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_peak_responses_%s.csv', exit_rule_suffix, timestamp));

share_common_peak_tex_path = fullfile(results_folder, ...
    sprintf('IRF_two_type_tauchen_%s_type_share_peak_responses_%s.tex', exit_rule_suffix, timestamp));

writetable(irf_out, csv_path);
writetable(rho_robustness_out, rho_robust_csv_path);
writetable(shock_size_robustness_out, shock_size_robust_csv_path);
writetable(share_robustness_out, share_robust_csv_path);
writetable(share_common_normalization_out, share_common_csv_path);
writetable(share_common_peak_table, share_common_peak_csv_path);
write_peak_response_latex_table(share_common_peak_tex_path, share_common_peak_table, shock_size, rho_pe);

save(mat_path, ...
    'years','Pe_index', ...
    'exit_rate_agg_pct','exit_rate_type_pct','exit_diff_pct', ...
    'mass_pre_total_path','mass_pre_type_path', ...
    'mass_active_total_index','mass_active_type_index', ...
    'entry_flow_total_path','entry_flow_type_path', ...
    'Y_index','Y_type_index', ...
    'E_index','E_type_index', ...
    'N_index','N_type_index', ...
    'impact_diff','peak_diff','peak_diff_idx','avg_diff_5y', ...
    'idx_low','idx_high', ...
    'shock_size','shock_size_robust_values','rho_pe','rho_robust_values','high_share_robust_values','chi_common','chi_type','zeta_A','smooth_window', ...
    'use_margin_exit_response','exit_margin_scale','exit_margin_gain', ...
    'baseline_exit_prob','exit_prob_path','vcheck_path','vcheck_base_policy');

save(rho_robust_mat_path, ...
    'rho_robust_values','rho_robust_irfs','rho_robustness_out', ...
    'shock_size','chi_common','chi_type','zeta_A','smooth_window', ...
    'use_margin_exit_response','exit_margin_scale','exit_margin_gain','idx_low','idx_high');

save(shock_size_robust_mat_path, ...
    'shock_size_robust_values','shock_size_robust_irfs','shock_size_robustness_out', ...
    'rho_pe','chi_common','chi_type','zeta_A','smooth_window', ...
    'use_margin_exit_response','exit_margin_scale','exit_margin_gain','idx_low','idx_high');

save(share_robust_mat_path, ...
    'high_share_robust_values','baseline_high_share','share_robust_irfs','share_robustness_out', ...
    'shock_size','rho_pe','chi_common','chi_type','zeta_A','smooth_window', ...
    'use_margin_exit_response','exit_margin_scale','exit_margin_gain','idx_low','idx_high');

save(share_common_mat_path, ...
    'high_share_robust_values','baseline_high_share','common_norm_denominators', ...
    'share_common_irfs','share_common_normalization_out','share_common_peak_table', ...
    'shock_size','rho_pe','idx_low','idx_high');

if exist('fig_symbolic','var') && isgraphics(fig_symbolic)
    try
        exportgraphics(fig_symbolic, symbolic_png_path, 'Resolution', 300);
        exportgraphics(fig_symbolic, symbolic_pdf_path, 'ContentType', 'vector');
    catch ME
        warning('Could not export symbolic IRF figure: %s', ME.message);
    end
end

if exist('fig_rho_robust','var') && isgraphics(fig_rho_robust)
    try
        exportgraphics(fig_rho_robust, rho_robust_png_path, 'Resolution', 300);
        exportgraphics(fig_rho_robust, rho_robust_pdf_path, 'ContentType', 'vector');
    catch ME
        warning('Could not export rho-robustness IRF figure: %s', ME.message);
    end
end

if exist('fig_shock_size_robust','var') && isgraphics(fig_shock_size_robust)
    try
        exportgraphics(fig_shock_size_robust, shock_size_robust_png_path, 'Resolution', 300);
        exportgraphics(fig_shock_size_robust, shock_size_robust_pdf_path, 'ContentType', 'vector');
    catch ME
        warning('Could not export shock-size robustness IRF figure: %s', ME.message);
    end
end

if exist('fig_share_robust','var') && isgraphics(fig_share_robust)
    try
        exportgraphics(fig_share_robust, share_robust_png_path, 'Resolution', 300);
        exportgraphics(fig_share_robust, share_robust_pdf_path, 'ContentType', 'vector');
    catch ME
        warning('Could not export type-share robustness IRF figure: %s', ME.message);
    end
end

if exist('fig_share_common','var') && isgraphics(fig_share_common)
    try
        exportgraphics(fig_share_common, share_common_png_path, 'Resolution', 300);
        exportgraphics(fig_share_common, share_common_pdf_path, 'ContentType', 'vector');
    catch ME
        warning('Could not export common-normalized type-share figure: %s', ME.message);
    end
end

function format_symbolic_axis(plot_end)
    ax = gca;
    box(ax, 'on');
    grid(ax, 'off');
    xlim(ax, [0 plot_end-1]);

    if isprop(ax, 'Toolbar')
        ax.Toolbar.Visible = 'off';
    end

    set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'LineWidth', 1.1, ...
        'TickDir', 'out');

    ax.XLabel.FontSize = 13;
    ax.XLabel.Units = 'normalized';
    ax.XLabel.Position(2) = -0.22;

    ax.YLabel.FontSize = 13;
    ax.Title.FontSize = 15;
end

function p_exit = marginal_exit_probability(vcheck, vcheck_base, hard_exit, scale, gain)
    p_base = smooth_exit_probability(vcheck_base, scale);
    p_shock = smooth_exit_probability(vcheck, scale);
    extra_exit = gain * max(p_shock - p_base, 0);
    p_exit = hard_exit + (1 - hard_exit) .* min(extra_exit, 1);
end

function p_exit = smooth_exit_probability(vcheck, scale)
    z = vcheck ./ scale;
    z = max(min(z, 40), -40);
    p_exit = 1 ./ (1 + exp(z));
end

function irf = simulate_rho_persistence_irf(rho_value, p)
    T = p.T;
    n_theta = p.n_theta;
    n_types = p.n_types;

    W_path = ones(T,1) * p.W;

    logPe_path = zeros(T,1);
    logPe_path(1) = log(1 + p.shock_size);

    for t = 2:T
        logPe_path(t) = rho_value * logPe_path(t-1);
    end

    Pe_path = p.Pe * exp(logPe_path);

    baseline_exit_prob = p.exit_policy;
    vcheck_base_policy = zeros(n_theta, n_types);

    for s = 1:n_types
        vcheck_base_policy(:,s) = -p.W*p.x + p.beta * (p.Q_incumbent * p.V_type(:,s));
    end

    V_path = zeros(n_theta, T, n_types);
    exit_policy_path = zeros(n_theta, T, n_types);
    exit_prob_path = zeros(n_theta, T, n_types);

    y_policy_path = zeros(n_theta, T, n_types);
    n_policy_path = zeros(n_theta, T, n_types);
    e_policy_path = zeros(n_theta, T, n_types);
    pi_policy_path = zeros(n_theta, T, n_types);

    V_tomorrow = p.V_type;

    for t = T:-1:1
        W_today = W_path(t);
        Pe_today = Pe_path(t);

        for s = 1:n_types
            gamma_s = p.gamma_type(s);
            A_s = p.A_type(s);
            shock_gap = max(Pe_today/p.Pe - 1, 0);

            A_today = A_s * (1 - p.zeta_A * shock_gap);
            A_today = max(A_today, 1e-6);

            eta_s = 1 - p.alpha - gamma_s;

            for i = 1:n_theta
                th = p.theta_grid(i);

                y_today = (A_today * th)^(1/eta_s) ...
                    * (p.alpha / W_today)^(p.alpha / eta_s) ...
                    * (gamma_s / Pe_today)^(gamma_s / eta_s);

                n_today = p.alpha * y_today / W_today;
                e_today = gamma_s * y_today / Pe_today;

                pi_today = y_today - W_today*n_today - Pe_today*e_today;

                y_policy_path(i,t,s) = y_today;
                n_policy_path(i,t,s) = n_today;
                e_policy_path(i,t,s) = e_today;
                pi_policy_path(i,t,s) = pi_today;
            end

            EV = p.beta * (p.Q_incumbent * V_tomorrow(:,s));
            gamma_rel = gamma_s / mean(p.gamma_type);
            x_eff = p.x * (1 + p.chi_common*shock_gap + p.chi_type*shock_gap*gamma_rel);

            vcheck = -W_today*x_eff + EV;
            exit_policy_path(:,t,s) = double(vcheck <= 0);

            if p.use_margin_exit_response
                exit_prob_path(:,t,s) = marginal_exit_probability( ...
                    vcheck, vcheck_base_policy(:,s), p.exit_policy(:,s), ...
                    p.exit_margin_scale, p.exit_margin_gain);
            else
                exit_prob_path(:,t,s) = exit_policy_path(:,t,s);
            end

            V_path(:,t,s) = pi_policy_path(:,t,s) + max(vcheck, 0);
        end

        V_tomorrow = squeeze(V_path(:,t,:));
    end

    m0_pre = p.M_total * p.mu;

    y0_policy = zeros(n_theta, n_types);
    n0_policy = zeros(n_theta, n_types);
    e0_policy = zeros(n_theta, n_types);
    pi0_policy = zeros(n_theta, n_types);

    for s = 1:n_types
        gamma_s = p.gamma_type(s);
        A_s = p.A_type(s);
        eta_s = 1 - p.alpha - gamma_s;

        for i = 1:n_theta
            th = p.theta_grid(i);

            y0_policy(i,s) = (A_s * th)^(1/eta_s) ...
                * (p.alpha / p.W)^(p.alpha / eta_s) ...
                * (gamma_s / p.Pe)^(gamma_s / eta_s);

            n0_policy(i,s) = p.alpha * y0_policy(i,s) / p.W;
            e0_policy(i,s) = gamma_s * y0_policy(i,s) / p.Pe;
            pi0_policy(i,s) = y0_policy(i,s) - p.W*n0_policy(i,s) - p.Pe*e0_policy(i,s);
        end
    end

    m0_operate = zeros(n_theta, n_types);

    for s = 1:n_types
        m0_operate(:,s) = m0_pre(:,s) .* (1 - baseline_exit_prob(:,s));
    end

    M0_type_active = zeros(1,n_types);
    Y0_type = zeros(1,n_types);
    E0_type = zeros(1,n_types);

    for s = 1:n_types
        M0_type_active(s) = sum(m0_operate(:,s));
        Y0_type(s) = sum(y0_policy(:,s) .* m0_operate(:,s));
        E0_type(s) = sum(e0_policy(:,s) .* m0_operate(:,s));
    end

    Y0_agg = sum(Y0_type);
    E0_agg = sum(E0_type);

    m_current = m0_pre;

    mass_pre_total_path = zeros(T+1,1);
    mass_pre_type_path = zeros(T+1,n_types);
    mass_active_type_path = zeros(T+1,n_types);

    exit_flow_total_path = zeros(T,1);
    exit_flow_type_path = zeros(T,n_types);
    entry_flow_type_path = zeros(T,n_types);

    exit_rate_type_path = zeros(T+1,n_types);

    Y_agg_path = zeros(T+1,1);
    E_agg_path = zeros(T+1,1);

    Y_type_path = zeros(T+1,n_types);
    E_type_path = zeros(T+1,n_types);

    mass_pre_total_path(1) = sum(m0_pre(:));

    for s = 1:n_types
        mass_pre_type_path(1,s) = sum(m0_pre(:,s));
    end

    mass_active_type_path(1,:) = M0_type_active;
    exit_rate_type_path(1,:) = p.exit_rate_type(:)';

    Y_agg_path(1) = Y0_agg;
    E_agg_path(1) = E0_agg;

    Y_type_path(1,:) = Y0_type;
    E_type_path(1,:) = E0_type;

    for t = 1:T
        mass_pre_total_path(t+1) = sum(m_current(:));

        for s = 1:n_types
            mass_pre_type_path(t+1,s) = sum(m_current(:,s));
        end

        m_operate = zeros(n_theta, n_types);
        m_next = zeros(n_theta, n_types);

        for s = 1:n_types
            exit_prob_t = exit_prob_path(:,t,s);
            operate = 1 - exit_prob_t;
            m_operate(:,s) = m_current(:,s) .* operate;

            Y_type_path(t+1,s) = sum(y_policy_path(:,t,s) .* m_operate(:,s));
            E_type_path(t+1,s) = sum(e_policy_path(:,t,s) .* m_operate(:,s));

            exit_flow_type_path(t,s) = sum(m_current(:,s) .* exit_prob_t);

            if mass_pre_type_path(t+1,s) > 0
                exit_rate_type_path(t+1,s) = exit_flow_type_path(t,s) / mass_pre_type_path(t+1,s);
            else
                exit_rate_type_path(t+1,s) = 0;
            end

            entry_flow_type_t = exit_flow_type_path(t,s);
            entry_flow_type_path(t,s) = entry_flow_type_t;

            m_next(:,s) = p.Q_incumbent' * m_operate(:,s) + entry_flow_type_t * p.Q_entrant;
        end

        for s = 1:n_types
            mass_active_type_path(t+1,s) = sum(m_operate(:,s));
        end

        Y_agg_path(t+1) = sum(Y_type_path(t+1,:));
        E_agg_path(t+1) = sum(E_type_path(t+1,:));

        exit_flow_total_path(t) = sum(exit_flow_type_path(t,:));
        m_next(m_next < 0) = 0;
        m_current = m_next;
    end

    exit_rate_agg_path = zeros(T+1,1);
    exit_rate_agg_path(1) = sum(mass_pre_type_path(1,:) .* p.exit_rate_type(:)') ...
        / mass_pre_total_path(1);

    valid_mass = mass_pre_total_path(2:end) > 0;
    exit_rate_agg_path(1+find(valid_mass)) = ...
        exit_flow_total_path(valid_mass) ./ mass_pre_total_path(1+find(valid_mass));

    years = 0:T;

    irf = struct();
    irf.rho_pe = rho_value;
    irf.shock_size = p.shock_size;
    type_shares = sum(p.mu,1) / sum(p.mu(:));
    irf.high_type_share = type_shares(p.idx_high);
    irf.years = years(:);
    irf.Y_ss = Y0_agg;
    irf.Y_low_ss = Y0_type(p.idx_low);
    irf.Y_high_ss = Y0_type(p.idx_high);
    irf.E_ss = E0_agg;
    irf.E_low_ss = E0_type(p.idx_low);
    irf.E_high_ss = E0_type(p.idx_high);
    irf.mass_active_ss = sum(M0_type_active);
    irf.mass_active_low_ss = M0_type_active(p.idx_low);
    irf.mass_active_high_ss = M0_type_active(p.idx_high);
    irf.Y_index = 100 * Y_agg_path / Y0_agg;
    irf.Y_low_index = 100 * Y_type_path(:,p.idx_low) / Y0_type(p.idx_low);
    irf.Y_high_index = 100 * Y_type_path(:,p.idx_high) / Y0_type(p.idx_high);
    irf.E_index = 100 * E_agg_path / E0_agg;
    irf.E_low_index = 100 * E_type_path(:,p.idx_low) / E0_type(p.idx_low);
    irf.E_high_index = 100 * E_type_path(:,p.idx_high) / E0_type(p.idx_high);
    irf.Y_diff_high_minus_low_indexpts = irf.Y_high_index - irf.Y_low_index;
    irf.E_diff_high_minus_low_indexpts = irf.E_high_index - irf.E_low_index;
    irf.exit_rate_agg_pct = 100 * exit_rate_agg_path;
    irf.exit_rate_low_pct = 100 * exit_rate_type_path(:,p.idx_low);
    irf.exit_rate_high_pct = 100 * exit_rate_type_path(:,p.idx_high);
    irf.exit_diff_high_minus_low_pctpts = irf.exit_rate_high_pct - irf.exit_rate_low_pct;
    irf.mass_active_index = 100 * sum(mass_active_type_path,2) / sum(M0_type_active);
    irf.mass_active_low_index = 100 * mass_active_type_path(:,p.idx_low) / M0_type_active(p.idx_low);
    irf.mass_active_high_index = 100 * mass_active_type_path(:,p.idx_high) / M0_type_active(p.idx_high);
    irf.mass_active_diff_high_minus_low_index = irf.mass_active_high_index - irf.mass_active_low_index;
end

function colors = rho_persistence_colors(n_colors)
    base_colors = [ ...
        0.00 0.00 0.00; ...
        0.00 0.20 0.70; ...
        0.70 0.05 0.10; ...
        0.20 0.50 0.20; ...
        0.55 0.25 0.60; ...
        0.85 0.45 0.05];

    if n_colors <= size(base_colors,1)
        colors = base_colors(1:n_colors,:);
    else
        colors = lines(n_colors);
        colors(1:size(base_colors,1),:) = base_colors;
    end
end

function add_rho_persistence_legend(fig_handle, rho_values, colors)
    add_value_legend(fig_handle, rho_value_labels(rho_values), colors);
end

function labels = rho_value_labels(rho_values)
    labels = cell(1,numel(rho_values));

    for ii = 1:numel(rho_values)
        labels{ii} = sprintf('$\\rho_{P_e} = %.2f$', rho_values(ii));
    end
end

function add_shock_size_legend(fig_handle, shock_values, colors)
    add_value_legend(fig_handle, shock_size_value_labels(shock_values), colors);
end

function labels = shock_size_value_labels(shock_values)
    labels = cell(1,numel(shock_values));

    for ii = 1:numel(shock_values)
        labels{ii} = sprintf('%.0f%% shock', 100*shock_values(ii));
    end
end

function add_high_share_legend(fig_handle, high_share_values, colors)
    add_value_legend(fig_handle, high_share_value_labels(high_share_values), colors);
end

function add_firm_type_legend(fig_handle)
    firm_type_colors = [ ...
        0.00 0.00 0.00; ...
        0.00 0.20 0.70; ...
        0.70 0.05 0.10];
    add_value_legend(fig_handle, {'Aggregate','Low-energy firms','High-energy firms'}, firm_type_colors);
end

function add_firm_type_column_headers(fig_handle)
    annotation(fig_handle, 'textbox', [0.115 0.925 0.33 0.028], ...
        'String', 'Low-energy firms', ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'FontWeight', 'bold');
    annotation(fig_handle, 'textbox', [0.595 0.925 0.33 0.028], ...
        'String', 'High-energy firms', ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'FontWeight', 'bold');
end

function labels = high_share_value_labels(high_share_values)
    labels = cell(1,numel(high_share_values));

    for ii = 1:numel(high_share_values)
        labels{ii} = sprintf('$s_H = %.0f\\%%$', 100*high_share_values(ii));
    end
end

function mu_reweighted = reweight_type_shares(mu_base, idx_low, idx_high, high_share)
    high_share = max(min(high_share, 0.95), 0.05);
    target_shares = zeros(1,size(mu_base,2));
    target_shares(idx_low) = 1 - high_share;
    target_shares(idx_high) = high_share;

    mu_reweighted = zeros(size(mu_base));
    total_mu = sum(mu_base(:));

    for s = 1:size(mu_base,2)
        col_sum = sum(mu_base(:,s));

        if col_sum > 0
            mu_reweighted(:,s) = total_mu * target_shares(s) * mu_base(:,s) / col_sum;
        end
    end
end

function add_value_legend(fig_handle, value_labels, colors)
    add_legend_group(fig_handle, '', value_labels, colors, [0.255 0.058 0.490 0.055]);
end

function add_dual_value_legend(fig_handle, left_title, left_labels, left_colors, right_title, right_labels, right_colors)
    add_legend_group(fig_handle, left_title, left_labels, left_colors, [0.040 0.058 0.455 0.055]);
    add_legend_group(fig_handle, right_title, right_labels, right_colors, [0.505 0.058 0.455 0.055]);
end

function add_legend_group(fig_handle, group_title, value_labels, colors, box_position)
    n_items = numel(value_labels);
    annotation(fig_handle, 'rectangle', box_position, ...
        'Color', [0 0 0], 'LineWidth', 1.0);

    title_height = 0.024;
    if ~isempty(group_title)
        annotation(fig_handle, 'textbox', [box_position(1)+0.008 box_position(2)+box_position(4)-0.031 box_position(3)-0.016 title_height], ...
            'String', group_title, ...
            'Interpreter', 'tex', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end

    x_left = box_position(1) + 0.016;
    x_right = box_position(1) + box_position(3) - 0.010;
    item_width = (x_right - x_left) / max(n_items,1);
    line_width = min(0.022, 0.17*item_width);
    legend_y = box_position(2) + 0.028;

    for ii = 1:n_items
        x0 = x_left + (ii-1)*item_width;
        annotation(fig_handle, 'line', [x0 x0+line_width], [legend_y legend_y], ...
            'Color', colors(ii,:), 'LineWidth', 1.8);
        label_interpreter = 'tex';
        if ischar(value_labels{ii}) && ~isempty(value_labels{ii}) && value_labels{ii}(1) == '$'
            label_interpreter = 'latex';
        end
        annotation(fig_handle, 'textbox', [x0+line_width+0.006 legend_y-0.017 max(0.110,item_width-line_width-0.008) 0.034], ...
            'String', value_labels{ii}, ...
            'Interpreter', label_interpreter, ...
            'EdgeColor', 'none', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 12, ...
            'VerticalAlignment', 'middle');
    end
end

function ax = plot_composition_scenario_panel(irf, aggregate_field, low_field, high_field, ...
        steady_state_values, y_label, panel_title, shock_plot_idx, tgrid_raw, tgrid_fine, plot_end, show_xlabel)
    ax = nexttile;
    hold(ax, 'on');

    fields = {aggregate_field, low_field, high_field};
    colors = [ ...
        0.00 0.00 0.00; ...
        0.00 0.20 0.70; ...
        0.70 0.05 0.10];

    for ii = 1:numel(fields)
        series = irf.(fields{ii});
        series_dev_fine = interp1(tgrid_raw, ...
            series(shock_plot_idx) - steady_state_values(ii), tgrid_fine, 'pchip');
        plot(ax, tgrid_fine, series_dev_fine, ...
            'Color', colors(ii,:), 'LineWidth', 1.8);
    end

    yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
    title(ax, panel_title, 'Interpreter','latex', ...
        'FontWeight','normal', 'FontSize', 15);
    ylabel(ax, y_label, 'FontName','Times New Roman', 'FontSize', 13);

    if show_xlabel
        xlabel(ax, 'Years', 'FontName','Times New Roman', 'FontSize', 13);
    end

    format_symbolic_axis(plot_end);
end

function plot_type_robustness_panel(irfs, colors, shock_plot_idx, tgrid_raw, tgrid_fine, field_name, steady_state_value, y_label, panel_title, plot_end, show_xlabel)
    ax = nexttile;
    hold(ax, 'on');

    for ii = 1:numel(irfs)
        series = irfs(ii).(field_name);
        series_dev_fine = interp1(tgrid_raw, series(shock_plot_idx) - steady_state_value, ...
            tgrid_fine, 'pchip');
        plot(ax, tgrid_fine, series_dev_fine, ...
            'Color', colors(ii,:), 'LineWidth', 1.8);
    end

    yline(ax, 0, 'k:', 'LineWidth', 1.1);
    title(ax, panel_title, 'Interpreter','latex', ...
        'FontWeight','normal', 'FontSize', 15);
    ylabel(ax, y_label, 'FontName','Times New Roman', 'FontSize', 13);

    if show_xlabel
        xlabel(ax, 'Years', 'FontName','Times New Roman', 'FontSize', 13);
    end

    format_symbolic_axis(plot_end);
end

function plot_diff_robustness_panel(irfs, colors, shock_plot_idx, tgrid_raw, tgrid_fine, field_name, y_label, panel_title, plot_end, show_xlabel)
    ax = nexttile;
    hold(ax, 'on');

    for ii = 1:numel(irfs)
        series = irfs(ii).(field_name);
        series_fine = interp1(tgrid_raw, series(shock_plot_idx), tgrid_fine, 'pchip');
        plot(ax, tgrid_fine, series_fine, ...
            'Color', colors(ii,:), 'LineWidth', 1.8);
    end

    yline(ax, 0, 'k:', 'LineWidth', 1.1);
    title(ax, panel_title, 'Interpreter','latex', ...
        'FontWeight','normal', 'FontSize', 15);
    ylabel(ax, y_label, 'FontName','Times New Roman', 'FontSize', 13);

    if show_xlabel
        xlabel(ax, 'Years', 'FontName','Times New Roman', 'FontSize', 13);
    end

    format_symbolic_axis(plot_end);
end

function peak_value = signed_peak_response(series)
    series = series(:);
    [~, peak_idx] = max(abs(series));
    peak_value = series(peak_idx);
end

function write_peak_response_latex_table(file_path, peak_table, shock_size, rho_pe)
    [file_id, error_message] = fopen(file_path, 'w');
    if file_id < 0
        error('Could not create LaTeX peak-response table: %s', error_message);
    end
    close_file = onCleanup(@() fclose(file_id));

    fprintf(file_id, '\\begin{table}[htbp]\n');
    fprintf(file_id, '\\centering\n');
    fprintf(file_id, '\\small\n');
    fprintf(file_id, '\\caption{Peak aggregate responses by high-energy firm share}\n');
    fprintf(file_id, '\\label{tab:type-share-peak-responses}\n');
    fprintf(file_id, '\\begin{minipage}{0.9\\textwidth}\n');
    fprintf(file_id, '\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrr@{}}\n');
    fprintf(file_id, '\\toprule\n');
    fprintf(file_id, 'High-energy share & Peak $Y$ & Peak $\\delta$ & Peak $M$ & Peak $E$ \\\\\n');
    fprintf(file_id, '\\midrule\n');

    for ii = 1:height(peak_table)
        fprintf(file_id, '%s & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
            peak_table.share_label{ii}, peak_table.peak_Y_pct(ii), ...
            peak_table.peak_delta_pct(ii), peak_table.peak_M_pct(ii), ...
            peak_table.peak_E_pct(ii));
    end

    fprintf(file_id, '\\bottomrule\n');
    fprintf(file_id, '\\end{tabular*}\n');
    fprintf(file_id, '\\par\\vspace{0.3em}\n');
    fprintf(file_id, '{\\footnotesize\\raggedright \\textit{Notes:} ');
    fprintf(file_id, 'The energy-price shock is %.0f\\%%, with $\\rho_{P_e}=%.2f$. ', ...
        100*shock_size, rho_pe);
    fprintf(file_id, ['Responses are percentages of the calibrated-composition steady state. ' ...
        'Each entry is the signed response with the largest absolute magnitude over the plotted horizon.\\par}\n']);
    fprintf(file_id, '\\end{minipage}\n');
    fprintf(file_id, '\\end{table}\n');

    clear close_file
end

%Notes for 29th Juni%Notes for 29th Juni
%Does it matter the periodicity?
%Is it correct the approach for a longer time period recovery with a higher
%persistence?
