% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% x=sthresh(x,gamma)
%
% This function performs soft threshholding on a vector x.  The
% result has
%
%  x(i)= {x(i)-gamma           x(i) > gamma
%        { 0                   -gamma <= x(i) <= gamma
%        {x(i)+gamma           x(i)< -gamma
%
% Note that this function works on the vector x "in place" for
% greater speed.
%
function x = sthresh(x, gamma)
%
% First, take gamma off the abs(x).
%
xs = abs(x) - gamma;
xs(xs < 0) = 0;
%
% Put back the sign.
%
x = xs .* sign(x);
