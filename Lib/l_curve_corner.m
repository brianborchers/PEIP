% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,reg_param)
% returns l curve corner estimated using a maximum curvature (kappa) estimation
% in log-log space
% rho is the misfit and eta is the model norm or seminorm
%
% INPUT
%   rho       - misfit
%   eta       - model norm or seminorm
%   reg_param - the regularization parameter
%
% OUTPUT
%   reg_corner  - the value of reg_param with maximum curvature
%   ireg_corner - the index of the value in reg_param with maximum curvature
%   kappa       - the curvature for each reg_param
%
function [reg_corner, ireg_corner, kappa] = l_curve_corner(rho, eta, reg_param)

%transform rho and eta into log-log space
x = log(rho);
y = log(eta);

% Triangular/circumscribed circle simple approximation to curvature
% (after Roger Stafford)

% the series of points used for the triangle/circle
x1 = x(1:end-2);
x2 = x(2:end-1);
x3 = x(3:end);
y1 = y(1:end-2);
y2 = y(2:end-1);
y3 = y(3:end);

% the side lengths for each triangle
a = sqrt((x3 - x2).^2+(y3 - y2).^2);
b = sqrt((x1 - x3).^2+(y1 - y3).^2);
c = sqrt((x2 - x1).^2+(y2 - y1).^2);

s = (a + b + c) / 2; %semi-perimeter

% the radius of each circle
R = (a .* b .* c) ./ (4 * sqrt((s .* (s - a) .* (s - b) .* (s - c))));

% The curvature for each estimate for each value which is
% the reciprocal of its circumscribed radius. Since there aren't circles for
% the end points they have no curvature
kappa = [0; 1 ./ R; 0];
[~, ireg_corner] = max(abs(kappa(2:end-1)));
reg_corner = reg_param(ireg_corner);
return;
