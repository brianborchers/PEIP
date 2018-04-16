% Example 9.2 forward problem
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% computes y=p(1)*exp(p(2)*x)+p(3)*x*exp(p(4)*x)
% 
% y=rawfunc(x,p)
%
function y=rawfunc(x,p)
y=p(1).*exp(p(2)*x)+p(3).*x.*exp(p(4)*x);
