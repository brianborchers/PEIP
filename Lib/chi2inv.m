% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=chi2inv(p,nu)
%
% Computes the inverse Chi^2 distribution corresponding to a given
% probability that a Chi^2 random variable with the given degrees
% of freedom is less than or equal to x.  Uses chi2cdf.m.
%
% Input Parameters:
%   p - probability that Chi^2 random variable is less than or
%       equal to x (scalar).
%   nu - degrees of freedom (scalar)
%
% Output Parameters:
%   x - corresponding value of x for given probability.
%
% Note that x and m must be scalars.
%
function x = chi2inv(p, nu)

% Special cases.
if (p >= 1.0)
    x = +Inf;
    return;
elseif (p == 0.0)
    x = 0;
    return
elseif (p < 0)
    x = -Inf;
    return;
end

% find a window with a cdf containing p
l = 0.0;
r = 1.0;
while (chi2cdf(r, nu) < p)
    l = r;
    r = r * 2;
end

% do a binary search until we have a sufficiently small interval around x
while (((r - l) / r) > 1.0e-5)
    m = (l + r) / 2;
    if (chi2cdf(m, nu) > p)
        r = m;
    else
        l = m;
    end
end
x = (l + r) / 2;
