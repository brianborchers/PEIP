% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
%  [a,u,b]=bracket(@f)
%
%  For a unimodal function f(alpha), finds a bracket [a,u,b] surrounding a
%  minimum of f(alpha), where f(a) >= f(u) and f(u) <= f(b).
%
function [a,u,b]=bracket(fun)
f0=feval(fun,0);
f1=feval(fun,1);
x=0.5;
if (f1 > f0) 
  while (feval(fun,x) > f0) 
    x=x/2;
  end
  a=0;
  b=2*x;
  u=x;
else
%
% The bracket will be bigger than 0,1
%
  x=2;
  while (feval(fun,x) < f1)
    x=2*x;
  end
  a=0;
  u=x/2;
  b=x;
end
