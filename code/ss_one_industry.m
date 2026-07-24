% Dynamic Macro Summer 2025
% Lab 6
% Partial equilibrium Germany design
% Unified one-industry baseline script with selectable grid method
clear; clc;

code_dir = fileparts(mfilename('fullpath'));
root = fileparts(code_dir);
addpath(code_dir, fullfile(code_dir, 'helpers'));
baselines_dir = fullfile(root, 'data', 'baselines');

% Select discretization path: 'tauchen' or 'rouwenhorst'
grid_method = 'tauchen';

% Parameters
beta        = 0.9500;
x           = 0.4600;   % fixed operating cost 
xE          = 35.17691728;   % entry cost 
rho_theta   = 0.9500;
sigma_theta = 0.1200;         
alpha       = 0.22;
gamma       = 0.016;
eta         = 1 - alpha - gamma;
L_supply    = 1;

% Iterate on W (1) or treat xE as a free parameter (0)?
W_iter_ind  = 0;   % =1 solve vE(W)=0, =0 recover xE from free entry

% Productivity grid
n_theta = 101;
switch lower(grid_method)
    case 'tauchen'
        [z, Q_incumbent] = tauchen(n_theta,0,rho_theta,sigma_theta,4);
    case 'rouwenhorst'
        [z, Q_incumbent] = rouwenhorst_ar1(n_theta,0,rho_theta,sigma_theta);
    otherwise
        error('Unknown grid_method. Use ''tauchen'' or ''rouwenhorst''.');
end
theta_grid = exp(z);

% Entrant distribution = stationary distribution
[V,D] = eig(Q_incumbent');
[~,idx] = min(abs(diag(D)-1));
pi_stat = real(V(:,idx));
pi_stat = pi_stat / sum(pi_stat);
Q_entrant = pi_stat;
Q_entrant = Q_entrant / sum(Q_entrant);

clear z V D idx pi_stat

% Wage iteration specifics
W_iter      = 0;
W_max_iter  = 1000;
W_diff      = 100;
W_conv_diff = 1e-6;
W_eps       = 0.01;

% Initial guess
W  = 1;
Pe = 1;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER AGGREGATE CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (W_iter < W_max_iter) && (W_diff > W_conv_diff)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATIC POLICIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_policy  = zeros(n_theta,1);
    e_policy  = zeros(n_theta,1);
    y_policy  = zeros(n_theta,1);
    pi_policy = zeros(n_theta,1);

    for theta_ind = 1:n_theta
        theta = theta_grid(theta_ind,1);

        y_policy(theta_ind,1) = theta^(1/eta) ...
                              * (alpha/W)^(alpha/eta) ...
                              * (gamma/Pe)^(gamma/eta);

        n_policy(theta_ind,1) = alpha * y_policy(theta_ind,1) / W;
        e_policy(theta_ind,1) = gamma * y_policy(theta_ind,1) / Pe;

        pi_policy(theta_ind,1) = y_policy(theta_ind,1) ...
                               - W * n_policy(theta_ind,1) ...
                               - Pe * e_policy(theta_ind,1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VFI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_diff      = 100;
    v_conv_diff = 1e-15;
    v_iter      = 0;
    v_max_iter  = 10000;

    v0          = zeros(n_theta,1);
    v1          = zeros(n_theta,1);
    vcheck      = zeros(n_theta,1);
    vhat        = zeros(n_theta,1);
    exit_policy = zeros(n_theta,1);

    while (v_diff > v_conv_diff) && (v_iter < v_max_iter)

        for theta_ind = 1:n_theta

            vp_exp = 0;
            for thetap_ind = 1:n_theta
                vp_exp = vp_exp + beta * Q_incumbent(theta_ind,thetap_ind) * v0(thetap_ind,1);
            end

            % Fixed operating cost in labor units
            vcheck(theta_ind,1) = -W * x + vp_exp;

            if vcheck(theta_ind,1) > 0
                vhat(theta_ind,1)        = vcheck(theta_ind,1);
                exit_policy(theta_ind,1) = 0;
            else
                vhat(theta_ind,1)        = 0;
                exit_policy(theta_ind,1) = 1;
            end

            v1(theta_ind,1) = pi_policy(theta_ind,1) + vhat(theta_ind,1);
        end

        v_diff = max(abs(v1 - v0));
        v_iter = v_iter + 1;
        v0     = v1;
    end

    if v_iter == v_max_iter
            warning (' VFI hit max iterations without full convergence');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VALUE TO ENTRY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vE = 0;
    for theta_ind = 1:n_theta
        vE = vE + beta * Q_entrant(theta_ind,1) * v0(theta_ind,1);
    end

    if W_iter_ind == 0
        % free parameter xE recovered below
        vE = vE;
    else
        % xE fixed in labor units
        vE = -W * xE + vE;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UPDATE W OR xE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if W_iter_ind == 1
        W_diff_signed = vE;
        W_diff        = abs(W_diff_signed);
        W_iter        = W_iter + 1;

        if W_diff > W_conv_diff
            W = W + W_eps * W_diff_signed;
        end
    else
        % recover xE in labor units
        xE     = vE / W;
        W_diff = 0;
    end

end

disp('W_diff, v_diff');
fprintf('W_diff = %.6e, v_diff = %.6e\n', W_diff, v_diff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROSS-SECTIONAL MEASURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_diff      = 100;
cs_iter      = 0;
cs_max_iter  = 20000;
cs_conv_diff = 1e-12;

mu = ones(n_theta,1) / n_theta;

while (cs_diff > cs_conv_diff) && (cs_iter < cs_max_iter)

    % 1. Surviving firms
    mu_survive = mu .* (1 - exit_policy);
    exit_rate_tmp = 1 - sum(mu_survive);

    % 2. Transition of surviving incumbents
    mu_next_inc = Q_incumbent' * mu_survive;

    % 3. Entry = exit in steady state
    s_e = exit_rate_tmp;
    mu_next = mu_next_inc + s_e * Q_entrant;
    mu_next = mu_next / sum(mu_next);

    cs_diff = max(abs(mu_next - mu));
    cs_iter = cs_iter + 1;
    mu      = mu_next;

end

if cs_iter == cs_max_iter
    warning('Stationary dist. iteration hit max iterations without full convergence.');
end

mu_survive   = mu .* (1 - exit_policy);
exit_rate    = max(0, 1 - sum(mu_survive));
mu_underline = Q_incumbent' * mu;

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_ave_temp  = zeros(n_theta,1);
y_ave_temp  = zeros(n_theta,1);
e_ave_temp  = zeros(n_theta,1);
pi_ave_temp = zeros(n_theta,1);

for theta_ind = 1:n_theta
    e_ave_temp(theta_ind,1)  = e_policy(theta_ind,1)  * mu_underline(theta_ind,1);
    n_ave_temp(theta_ind,1)  = n_policy(theta_ind,1)  * mu_underline(theta_ind,1);
    y_ave_temp(theta_ind,1)  = y_policy(theta_ind,1)  * mu_underline(theta_ind,1);
    pi_ave_temp(theta_ind,1) = pi_policy(theta_ind,1) * mu_underline(theta_ind,1);
end

e_ave  = sum(e_ave_temp);
n_ave  = sum(n_ave_temp);
y_ave  = sum(y_ave_temp);
pi_ave = sum(pi_ave_temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND MASS OF ACTIVE FIRMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Labor demand includes variable labor + fixed operating labor + entry labor
M          = L_supply / (n_ave + (1 - exit_rate)*x + exit_rate*xE);
M_entrants = M * exit_rate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGGREGATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_agg  = M * e_ave;
n_agg  = M * n_ave;
y_agg  = M * y_ave;
pi_agg = M * pi_ave;

% Fixed costs in goods value = wage times labor units
fc_agg = W * ( (M - M_entrants) * x + M_entrants * xE );

pi_var_agg        = pi_agg;
pi_net_agg        = pi_agg - fc_agg;
profit_margin_net = pi_net_agg / y_agg;
profit_margin_var = pi_var_agg / y_agg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Household consumption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy revenues are rebated lump-sum to households
C = Pe * e_agg + W * L_supply + pi_agg - fc_agg;

% Goods market consistency check
if abs(C - y_agg) > 1e-6
    disp('Goods market identity violation!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE BASELINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_initial  = v0;
mu_initial = mu;

baseline_file = fullfile(baselines_dir, ...
    sprintf('German_Baseline_OneType_%s.mat', upper(grid_method)));
save(baseline_file, ...
    'V_initial', 'mu_initial', ...
    'alpha', 'gamma', 'x', 'xE', 'beta', ...
    'theta_grid', 'n_theta', 'Q_incumbent', 'Q_entrant', 'eta', ...
    'exit_rate', 'M', 'M_entrants', ...
    'e_agg', 'L_supply', 'W', 'Pe', 'y_agg', 'n_agg', 'pi_agg', ...
    'grid_method', 'fc_agg');

fprintf('Saved one-industry baseline: %s\n', baseline_file);
