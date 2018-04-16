% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=phiinv(z)
%
% Calculates the inverse normal distribution from the value of the
% integral
%
%       z=int((1/sqrt(2*pi))*exp(-t^2/2),t=-infinity..x)
%
% and returns the value of x corresponding to the integral value.
%
% Input Parameters:
%   z - value of integral (scalar)
%
% Output Parameters:
%   x - endpoint value of integration (scalar)
function x = phiinv(z)

if (z >= 0.5),
    x = sqrt(2) * erfinv((z - 0.5)/.5);
else
    x = -phiinv(1-z);
end
