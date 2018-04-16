% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% l=loglikelihood(m)
%
function l=loglikelihood(m)
%
% global variables.
%
global x;
global y;
global sigma;
%
% Compute the standardized residuals.
%
fvec=(y-fun(m,x))./sigma;
%
% The log likelihood is (-1/2)*sum(fvec(i)^2,i=1..n);
%
l=(-1/2)*sum(fvec.^2);
