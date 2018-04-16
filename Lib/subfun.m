% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% f=subfun(alpha)
%
%  This computes 'FUN' at X+alpha*P.  Notice that FUN, X, and P are passed
%  in through global variables.  
%
function f=subfun(alpha)
global FUN;
global X;
global P;

f=feval(FUN,X+alpha*P);
