% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%   m=simmvn(mean,cov)
%
% Generates a random vector realization for  a specified multivariate normal distribution.
%
function m = simmvn(mean, cov)
n = length(mean);
R = chol(cov);
m = R' * randn(n, 1) + mean;
