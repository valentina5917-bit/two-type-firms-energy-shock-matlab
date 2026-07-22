function out = SS_moments(x, p)

% Steady-state moments for calibration of x and xE
% p.n_theta;% x and xE are fixed costs in labor units

n_theta = p.n_theta;

% Static firm policies
n_policy = zeros(n_theta,1);
y_policy = zeros(n_theta,1);

for i = 1:n_theta
    th = p.theta_grid(i);

    % Static profit maximization
    y = th^(1/p.eta) ...
        * (p.alpha/p.W)^(p.alpha/p.eta) ...
        * (p.gamma/p.Pe)^(p.gamma/p.eta);

    n_policy(i) = p.alpha * y / p.W;
    y_policy(i) = y;
end

e_policy  = p.gamma * y_policy / p.Pe;
pi_policy = y_policy - p.W*n_policy - p.Pe*e_policy;

% Incumbent value function iteration
V0 = zeros(n_theta,1);
exit_policy = false(n_theta,1);

for it = 1:20000
    EV = p.beta * (p.Q_incumbent * V0);

    % continuation value net of fixed operating cost in wage units
    vcheck = EV - p.W * x;

    V1 = pi_policy + max(vcheck, 0);
    exit_policy = (vcheck <= 0);

    if max(abs(V1 - V0)) < 1e-10
        break;
    end

    V0 = V1;
end

% Free-entry implied xE (in labor units)
% 0 = -W*xE + beta * E[V]
xE_FE = (p.beta * (p.Q_entrant' * V0)) / p.W;

% Stationary distribution with entry and exit
mu = ones(n_theta,1) / n_theta;

for it = 1:50000
    mu_survive = mu .* (1 - exit_policy);
    exit_rate = 1 - sum(mu_survive);

    mu_next = p.Q_incumbent' * mu_survive + exit_rate * p.Q_entrant;
    mu_next = mu_next / sum(mu_next);

    if max(abs(mu_next - mu)) < 1e-12
        mu = mu_next;
        break;
    end

    mu = mu_next;
end

exit_rate = max(exit_rate, 0);

% Outputs
out.exit        = exit_rate;   % yearly model exit probability
out.xE          = xE_FE;       % entry cost in labor units
out.V           = V0;          % diagnostic
out.exit_policy = exit_policy;

end


