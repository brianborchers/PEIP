% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% computes the forward problem y=p(1)*exp(p(2)*x)+p(3)*x*exp(p(4)*x)
% 
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
function y=fun(m,x)
y=m(1).*exp(m(2)*x)+m(3).*x.*exp(m(4)*x);
