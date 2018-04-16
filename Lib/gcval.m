% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% Evaluate the GCV function "gcv_function" at npoints points
%
% [reg_min,G,alpha] = gcval(U,s,b,npoints)
% s = [lambda,mu] evaluted from gsvd
%
% Returns the GCV function
% g = || Gm_(alpha,L) - d ||^2 / (Tr(I - GG#)^2
%
% INPUT
%   U       - the U matrix from gsvd(G, L)
%   s       - [diag(C) diag(S)] which are the lambdas and mus from the gsvd
%   b       - the data to try and match
%   npoints - the number of alpha to estimate
%
% OUTPUT
%   reg_min - the alpha with the minimal g (scalar)
%   g       - || Gm_(alpha,L) - d ||^2 / (Tr(I - GG#)^2
%   alpha   - the alpha for the corresponding g
function [reg_min, g, alpha] = gcval(U, s, b, npoints)

% Smallest regularization parameter.
smin_ratio = 16 * eps;

% get matrix sizes
[m, ~] = size(U);
[p, ~] = size(s);

% project the data into a more useful space
beta = U' * b;;

% sort s in the opposite order and divide the first element in each row by
% the second to evaluate gamma
gamma = s(p:-1:1, 1) ./ s(p:-1:1, 2);
beta = beta(m:-1:1);

% Vector of regularization parameters.
alpha = zeros(npoints, 1);
%initiatize g
g = zeros(npoints, 1);
gamma2 = gamma.^2;

% the last alpha is the larger of the smallest ratio, or the largest ratio
% times a very small number
alpha(npoints) = max([gamma(p), gamma(1) * smin_ratio]);
r = (gamma(1) / alpha(npoints))^(1 / (npoints - 1));
% logarithmically distribute the alphas to be used
for i = npoints - 1:-1:1
    alpha(i) = r * alpha(i+1);
end

% Evaluate GCV function values.
for i = 1:npoints
    g(i) = gcv_function(alpha(i), gamma2, beta);
end

% find the minimal value of g and save it into a output variable
[~, mingi] = min(g);
reg_min = alpha(mingi);
