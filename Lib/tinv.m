% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% t=tinv(p,nu)
%
% Computes the inverse t distribution value corresponding to the
% given probability and degrees of freedom.
%
% Input Parameters:
%   p - probability that a t random variable is less than or equal to x
%   (not vectorized)
%   nu - degrees of freedom
%
% Output Parameters:
%   t - value of t corresponding to the probability.
function t = tinv(p, nu)
%
% Special cases.
%
if (p >= 1.0),
    t = +Inf;
    return;
end
if (p <= 0.0),
    t = -Inf;
    return;
end
%
% First, figure out whether t>1/2 or t<1/2.
%
if (p < 0.5),
    t = -tinv(1-p, nu);
    return;
else
    l = 0.0;
    r = 1.0;
    while (tcdf(r, nu) < p)
        l = r;
        r = r * 2;
    end
    %
    % Now, we've got a bracket around t.
    %
    while (((r - l) / r) > 1.0e-5)
        m = (l + r) / 2;
        if (tcdf(m, nu) > p)
            r = m;
        else
            l = m;
        end
    end
end
t = (l + r) / 2;
