% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% return l curve parematers for Tikhonov Regularization
%
% Routine originally inspired by Per Hansen's l-curve
% program (http://www2.imm.dtu.dk/~pch/Regutools/)
%
%
% [rho,eta,reg_param] = l_curve_tikh(U,s,d,npoints,[alpha_min, alpha_max]) )
%
% INPUT
%   U       - matrix of data space basis vectors from the svd
%   s       - vector of singular values
%   d       - the data vector
%   npoints - the number of logarithmically spaced regularization parameters
%   [alpha_min, alpha_max] if specified, constrain the logrithmically spaced
%   regularization parameter range, otherwise an attempt is made to estimate
%   them from the range of singular values
%
% OUTPUT
%   eta       - the solution norm ||m|| or seminorm ||Lm||
%   rho       - the residual norm ||G m - d||
%   reg_param - corresponding regularization parameters
%

function [rho, eta, reg_param] = l_curve_tikh_svd(U, s, d, npoints, varargin)

% Initialization.
[m, n] = size(U);
[p] = length(s);

% compute the projection, and residual error introduced by the projection
d_proj = U' * d;
dr = norm(d)^2 - norm(d_proj)^2;

%data projections
d_proj = d_proj(1:p);

%scale series terms by singular values
d_proj_scale = d_proj ./ s;

% initialize storage space
eta = zeros(npoints, 1);
rho = eta;
reg_param = eta;
s2 = s.^2;

if size(varargin, 2) == 0
    % set the smallest regularization parameter that will be used
    smin_ratio = 16 * eps;
    reg_param(npoints) = max([s(p), s(1) * smin_ratio]);
    
    % ratio so that reg_param(1) will be s(1)
    ratio = (s(1) / reg_param(npoints))^(1 / (npoints - 1));
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

% determine the fit for each parameter
for i = 1:npoints
    %GSVD filter factors
    f = s2 ./ (s2 + reg_param(i)^2);
    eta(i) = norm(f.*d_proj_scale);
    rho(i) = norm((1 - f).*d_proj);
end

% if we couldn't match the data exactly add the projection induced misfit
if (m > n && dr > 0)
    rho = sqrt(rho.^2+dr);
end
