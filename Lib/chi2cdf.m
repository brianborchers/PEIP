% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% p=chi2cdf(x,m)
%
% Computes the Chi^2 CDF, using a transformation to N(0,1) on page
% 333 of Thistead, Elements of Statistical Computing.
%
% Input Parameters:
%   x - end value of chi^2 pdf to integrate to. (scalar)
%   m - degrees of freedom (scalar)
%
% Output Parameters:
%   p - probability that Chi^2 random variable is less than or
%       equal to x (scalar).
%
% Note that x and m must be scalars.
%
function p = chi2cdf(x, m)

if (x == (m - 1)),
    p = 0.5;
else
    z = (x - m + 2 / 3 - 0.08 / m) * sqrt((m - 1)*log((m - 1)/x)+x-(m - 1)) / abs(x-m+1);
    p = phi(z);
end
