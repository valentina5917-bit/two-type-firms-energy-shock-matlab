
function [z, P] = rouwenhorst_ar1(n, mu, rho, sigma_eps)
%-------------------------------------------------------%
%Rouwenhorst discretization for:
% z' = mu + rho*z + eps,   eps ~ N(0, sigma_eps^2)
%
%Inputs:
% n         = number of grid points
% mu        = constant term in AR(1)
% rho       = persistence
% sigma_eps = std. dev. of innovation
%
% Outputs:
% z = grid
% P = transition matrix
%--------------------------------------------------------%

    if n < 2
        error('n must be at least 2');
    end

    p = (1 + rho) / 2;
    q = p;

    %Base case
    P = [p, 1-p;
         1-q, q];

    %Recursive form
    for m = 3:n
        P_old = P;
        P = p     * [P_old, zeros(m-1,1); zeros(1,m-1), 0] ...
          + (1-p) * [zeros(m-1,1), P_old; 0, zeros(1,m-1)] ...
          + (1-q) * [zeros(1,m-1), 0; P_old, zeros(m-1,1)] ...
          + q     * [0, zeros(1,m-1); zeros(m-1,1), P_old];

        P(2:m-1, :) = P(2:m-1, :) / 2;
    end

    %Unconditional std of AR(1)
    sigma_z = sigma_eps / sqrt(1 - rho^2);

    %Rouwenhorst grid width
    psi = sqrt(n - 1) * sigma_z;

    %Unconditional mean
    zbar = mu / (1 - rho);

    z = linspace(-psi, psi, n)' + zbar;

    %Numerical cleanup
    P(P < 0) = 0;
    P = P ./ sum(P, 2);
end
