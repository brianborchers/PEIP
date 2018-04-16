% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% return l curve parematers for truncated gsvd tegularization
% Routine originally inspired by Per Hansen's l-curve
% program (http://www2.imm.dtu.dk/~pch/Regutools/)
%
%
% [rho,eta,reg_param] = l_curve_tgsvd(U,sm,d)
%
% INPUTS
%   gsvd U and sm [vector of lambda, mu], where mu is sorted in ascending order
%   d       - the data vector
%
% OUTPUT
%   eta       - the solution seminorm ||Lm||
%   rho       - the residual norm ||Gm-d||
%   reg_param - corresponding regularization parameters
%
function [rho, eta, reg_param] = l_curve_tsvd(U, sm, d)

% Initialization.
[m, n] = size(U);
[p, ps] = size(sm);

% compute the projection, and residual error introduced by the projection
d_proj = U' * d;
dr = norm(d)^2 - norm(d_proj)^2;

if (ps == 1)
    %SVD case
    s = sm;
    d_proj = d_proj(1:p);
else
    %GSVD case
    s = sm(p:-1:1, 1) ./ sm(p:-1:1, 2);
    d_proj = d_proj(p:-1:1);
end

%scale series terms by singular values
d_proj_scale = d_proj(1:p) ./ s;

% initialize storage space
eta = zeros(p, 1);
rho = eta;

% populate the solution (semi)norm
eta(1) = abs(d_proj_scale(1))^2;

% add to the solution (semi)norm for each singular value used
for k = 2:p
    eta(k) = eta(k-1) + abs(d_proj_scale(k))^2;
end

% make eta actually store the model (semi)norm
eta = sqrt(eta);

% populate the least noisy fit based on projection misfit
if (m > n)
    if (dr > 0)
        rho(p) = dr;
    else
        rho(p) = eps^2;
    end
else
    rho(p) = eps^2;
end

% add to the misfit for each singular value removed
for k = p - 1:-1:1
    rho(k) = rho(k+1) + abs(d_proj(k+1))^2;
end

% make rho actually store the misfit
rho = sqrt(rho);

reg_param = [1:p]';
