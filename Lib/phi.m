% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% z=phi(x)
%
% Calculates the normal distribution and returns the value of the
% integral
%
%       z=int((1/sqrt(2*pi))*exp(-t^2/2),t=-infinity..x)
%
% Input Parameters:
%   x - endpoint of integration (scalar)
%
% Output Parameters:
%   z - value of integral
function z = phi(x)

if (x >= 0)
    z = .5 + .5 * erf(x/sqrt(2));
else
    z = 1 - phi(-x);
end
