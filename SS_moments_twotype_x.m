function out = SS_moments_twotype_x(x, p)

n_theta = p.n_theta;
n_types = length(p.gamma_type);

alpha       = p.alpha;
beta        = p.beta;
W           = p.W;
Pe          = p.Pe;
gamma_type  = p.gamma_type;
lambda_type = p.lambda_type;

theta_grid  = p.theta_grid;
Q_incumbent = p.Q_incumbent;
Q_entrant   = p.Q_entrant;
L_supply    = p.L_supply;

%% 1. Static policies + VFI by type

eta_type    = zeros(1,n_types);
y_policy    = zeros(n_theta,n_types);
n_policy    = zeros(n_theta,n_types);
e_policy    = zeros(n_theta,n_types);
pi_policy   = zeros(n_theta,n_types);
V_type      = zeros(n_theta,n_types);
exit_policy = zeros(n_theta,n_types);

tol_vfi   = 1e-10;
maxit_vfi = 20000;

for s = 1:n_types

    gamma_s = gamma_type(s);
    eta_s   = 1 - alpha - gamma_s;
    eta_type(s) = eta_s;

    for i = 1:n_theta
        th = theta_grid(i);

        y_policy(i,s) = th^(1/eta_s) ...
                      * (alpha/W)^(alpha/eta_s) ...
                      * (gamma_s/Pe)^(gamma_s/eta_s);

        n_policy(i,s) = alpha * y_policy(i,s) / W;
        e_policy(i,s) = gamma_s * y_policy(i,s) / Pe;

        pi_policy(i,s) = y_policy(i,s) ...
                       - W * n_policy(i,s) ...
                       - Pe * e_policy(i,s);
    end

    V0 = zeros(n_theta,1);

    for it = 1:maxit_vfi
        EV = beta * (Q_incumbent * V0);
        vcheck = EV - W * x;

        V1 = pi_policy(:,s) + max(vcheck, 0);
        exit_policy(:,s) = (vcheck <= 0);

        if max(abs(V1 - V0)) < tol_vfi
            break;
        end

        V0 = V1;
    end

    V_type(:,s) = V1;
end

%% Free entry -> implied xE

vE_type = zeros(n_types,1);
for s = 1:n_types
    vE_type(s) = beta * (Q_entrant' * V_type(:,s));
end

vE_total = lambda_type * vE_type;
xE_implied = vE_total / W;

%% Stationary distribution across types

mu = zeros(n_theta,n_types);
for s = 1:n_types
    mu(:,s) = lambda_type(s) * ones(n_theta,1) / n_theta;
end
mu = mu / sum(mu(:));

tol_mu   = 1e-12;
maxit_mu = 50000;

for it_mu = 1:maxit_mu

    mu_next = zeros(n_theta,n_types);
    exit_mass_type = zeros(n_types,1);

    for s = 1:n_types
        mu_survive_s = mu(:,s) .* (1 - exit_policy(:,s));
        exit_mass_type(s) = sum(mu(:,s)) - sum(mu_survive_s);

        mu_next(:,s) = Q_incumbent' * mu_survive_s;
    end

    total_exit_mass = sum(exit_mass_type);

    for s = 1:n_types
        mu_next(:,s) = mu_next(:,s) + total_exit_mass * lambda_type(s) * Q_entrant;
    end

    mu_next(mu_next < 0) = 0;
    mu_next = mu_next / sum(mu_next(:));

    if max(abs(mu_next(:) - mu(:))) < tol_mu
        mu = mu_next;
        break;
    end

    mu = mu_next;
end

%% Type shares and exit rates

share_type = [sum(mu(:,1)); sum(mu(:,2))];

mu_survive = zeros(n_theta,n_types);
exit_rate_type = zeros(n_types,1);

for s = 1:n_types
    mu_survive(:,s) = mu(:,s) .* (1 - exit_policy(:,s));

    if sum(mu(:,s)) > 0
        exit_rate_type(s) = 1 - sum(mu_survive(:,s)) / sum(mu(:,s));
    else
        exit_rate_type(s) = 0;
    end
end

exit_rate_agg = share_type' * exit_rate_type;
exit_rate_agg = max(exit_rate_agg, 0);

%% Pre-exit / production distributon

mu_underline = zeros(n_theta,n_types);
for s = 1:n_types
    mu_underline(:,s) = Q_incumbent' * mu(:,s);
end

%% Type averages

y_ave_type  = zeros(n_types,1);
n_ave_type  = zeros(n_types,1);
e_ave_type  = zeros(n_types,1);
pi_ave_type = zeros(n_types,1);

for s = 1:n_types
    y_ave_type(s)  = sum(y_policy(:,s)  .* mu_underline(:,s));
    n_ave_type(s)  = sum(n_policy(:,s)  .* mu_underline(:,s));
    e_ave_type(s)  = sum(e_policy(:,s)  .* mu_underline(:,s));
    pi_ave_type(s) = sum(pi_policy(:,s) .* mu_underline(:,s));
end

y_ave_total  = sum(y_ave_type);
n_ave_total  = sum(n_ave_type);
e_ave_total  = sum(e_ave_type);
pi_ave_total = sum(pi_ave_type);

%% Mass of firms

M_total = L_supply / (n_ave_total + (1 - exit_rate_agg) * x + exit_rate_agg * xE_implied);
M_type  = M_total * share_type;

M_entrants_total = M_total * exit_rate_agg;
M_entrants_type  = M_entrants_total * lambda_type';

%% Aggregates 

y_agg = M_total * y_ave_total;
n_agg = M_total * n_ave_total;
e_agg = M_total * e_ave_total;
pi_agg = M_total * pi_ave_total;

Y_type = M_type .* y_ave_type;
N_type = M_type .* n_ave_type;
E_type = M_type .* e_ave_type;

fc_agg = W * ( (M_total - M_entrants_total) * x + M_entrants_total * xE_implied );

output_share_type = Y_type / sum(Y_type);
labor_share_type  = N_type / sum(N_type);
energy_share_type = E_type / sum(E_type);

%% Output
out.x                = x;
out.xE_implied       = xE_implied;

out.V_type           = V_type;
out.exit_policy      = exit_policy;
out.mu               = mu;
out.mu_underline     = mu_underline;

out.exit_rate_type   = exit_rate_type;
out.exit_rate_agg    = exit_rate_agg;

out.M_total          = M_total;
out.M_type           = M_type;
out.M_entrants_total = M_entrants_total;
out.M_entrants_type  = M_entrants_type;

out.y_ave_type       = y_ave_type;
out.n_ave_type       = n_ave_type;
out.e_ave_type       = e_ave_type;
out.pi_ave_type      = pi_ave_type;

out.y_agg            = y_agg;
out.n_agg            = n_agg;
out.e_agg            = e_agg;
out.pi_agg           = pi_agg;
out.fc_agg           = fc_agg;

out.output_share_type = output_share_type;
out.labor_share_type  = labor_share_type;
out.energy_share_type = energy_share_type;
end
