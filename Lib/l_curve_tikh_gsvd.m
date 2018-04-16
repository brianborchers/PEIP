% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% return l curve parameters and models for truncated gsvd regularization
%
% [rho,eta,reg_param,m] = l_curve_tikh2(U,d,X,Lam,Mu,G,L,npoints,[alpha_min, alpha_max])
%
% U, X, Lam, Mu are obtained from the gsvd of the system matrix G and the
% corresponding roughening matrix L
% d is the data vector for the problem G*m=d
% npoints is the nubmer of L-curve points returned
% [alpha_min, alpha_max] if specified, constrain the logrithmically spaced
% regularization parameter range, otherwise an attempt is made to estimate
% them from the range of generalized singular values
%
% OUTPUT
%   eta       - the solution seminorm ||Lm||
%   rho       - the residual norm ||G m - d||
%   reg_param - corresponding regularization parameters
%   m         - corresponding suite of models for truncated GSVD
%
%
function [rho, eta, reg_param, m] = l_curve_tikh_gsvd(U, d, X, Lam, Mu, G, L, npoints, varargin)

% Initialization.
[M, N] = size(G);
p = rank(L);
Y = inv(X)';
lambda = sqrt(diag(Lam'*Lam));
mu = sqrt(diag(Mu'*Mu));
dproj = U' * d;

gamma = lambda ./ mu;

m = zeros(N, npoints);
rho = zeros(npoints, 1);
eta = zeros(npoints, 1);

if size(varargin, 2) == 0
    smin_ratio = 16 * eps;
    
    % Minimum regularization parameter
    if (M <= N)
        reg_param(npoints) = max([gamma(p), gamma(N-M+1) * smin_ratio]);
        
        % ratio so that reg_param(1) will be s(1)
        ratio = (gamma(N-M+1) / reg_param(npoints))^(1 / (npoints - 1));
    else
        reg_param(npoints) = gamma(p);
        ratio = (gamma(N-M+1) / reg_param(npoints))^(1 / (npoints - 1));
    end
end

if size(varargin, 2) == 2
    alpharange = cell2mat(varargin);
    reg_param(npoints) = alpharange(2);
    ratio = (alpharange(1) / alpharange(2))^(1 / (npoints - 1));
end

% calculate all the regularization parameters
for i = npoints - 1:-1:1
    reg_param(i) = ratio * reg_param(i+1);
end

if M > N
    k = 0;
else
    k = N - M;
end
Y = inv(X)';
gamma = lambda ./ mu;
ng = length(gamma);

%solve for each solution
for i = 1:npoints
    %evaluate series filter coefficients for this regularization parameter
    f = zeros(ng, 1);
    for j = 1:ng
        f(j) = gamma(j)^2 / (gamma(j)^2 + reg_param(i)^2);
        if isnan(gamma(j)) || isinf(gamma(j))
            f(j) = 1;
        end
        if (lambda(j) == 0 && mu(j) == 0)
            f(j) = 0;
        end
        
        %build the solution
        
        m(:, i) = m(:, i) + f(j) * (U(:, j+k)' * d / lambda(j)) * Y(:, j);
    end
    
    rho(i) = norm(G*m(:, i)-d);
    eta(i) = norm(L*m(:, i));
end

