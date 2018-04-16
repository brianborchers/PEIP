% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% g=grad(u0)
%
function g=grad(u0)
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
% First, compute the function value and make sure that Fvec is set.
%
f=objfun(u0);
%
% Now, use the adjoint equation to find the gradient.  We want to
% take 
%
v=zeros(n,1);
v(xpoints)=2*Fvec;
g=adjoint(m,n,deltat,deltax,D,v);
%
% Add the gradient of the regularization term.
%
g=g+2*alpha*L'*(L*u0);
