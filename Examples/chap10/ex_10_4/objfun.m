% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% ssq=objfun(u0);
%
% Computes the sum of squares of the differences between u(x,T) and
% the data points.
% 
function ssq=objfun(u0);
%
% Global variables.
%
global m;
global n;
global deltax;
global x;
global deltat;
global xpoints;
global d;
global D;
global Fvec;
global alpha;
global L;
%
% Run the model forward.
%
ufinal=forward(m,n,deltat,deltax,D,u0);
%
% Compute the differences
%
Fvec=ufinal(xpoints)-d;
%
% Compute the sum of squares and add the regularization term.
%
ssq=norm(Fvec)^2+alpha*norm(L*u0)^2;
