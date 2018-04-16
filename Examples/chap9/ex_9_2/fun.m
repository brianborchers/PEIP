% Example 9.2. Computes the difference between the model described by p and the collected data
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
%  f=fun(p)
%
function f=fun(p)
% global variables, the x and y points and their sigmas
global x;
global y;
global sigma;

f=zeros(length(y),1);
f=(rawfunc(x,p)-y)./sigma;
