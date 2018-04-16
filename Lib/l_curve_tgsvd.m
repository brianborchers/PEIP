% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% return l curve parematers and models for truncated gsvd regularization
%
%
% [rho,eta,reg_param,m] = l_curve_tgsvd(U,d,X,Lam,G,L)
%
% OUTPUT
%   eta       - the solution seminorm ||Lm||
%   rho       - the residual norm ||G m - d||
%   reg_param - corresponding regularization parameters
%   m         - corresponding suite of models for truncated GSVD
%
%
function [rho, eta, reg_param, m] = l_curve_tgsvd(U, d, X, Lam, G, L)

% Initialization.
[M, N] = size(G);
p = rank(L);
reg_param = (1:p + 1)';
Y = inv(X)';
lambda = sqrt(diag(Lam'*Lam));
dproj = zeros(1, M);
for i = 1:M
    dproj(i) = U(:, i)' * d;
end

m = zeros(N, M);
%build the solutions
m(:, 1) = (dproj(end) / lambda(end)) * Y(:, end);
for q = 2:N
    m(:, q) = m(:, q-1) + (dproj(end-q+1) / lambda(end-q+1)) * Y(:, end-q+1);
end

% initialize output variables
eta = zeros(M, 1);
rho = eta;

% calculate the solution misfit and seminorm for each generalized singular value used
for i = 1:M
    rho(i) = norm(G*m(:, i)-d);
    eta(i) = norm(L*m(:, i));
end

