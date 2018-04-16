% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% u=adjoint(m,n,deltat,deltax,D,uT)
%
% Integrate the heat equation backward in time using the Crank-Nicolson
% implicit Euler method.
%
% The system of equations for each step in time is of the form
%
%    A*u^(n+1)=B*u^(n)
%
% where A and B are tridiagonal matrices.  
%
function [u,A,B]=adjoint(m,n,deltat,deltax,D,u0,A,B)
if (nargin < 8)
%
% Setup A.
%
e=ones(n,1);
A=spdiags([(-D*deltat/(2*deltax^2))*e (1+2*D*deltat/(2*deltax^2))*e ...
              (-D*deltat/(2*deltax^2))*e],-1:1,n,n);
%
% Setup B.
%
B=spdiags([((D*deltat)/(2*deltax^2))*e +(1-2*D*deltat/(2*deltax^2))*e ...
              (D*deltat/(2*deltax^2))*e],-1:1,n,n);
end
%
% Do the time steps.
%
u=u0;
for i=1:m
  u=(B')*(A'\u);
end
%
% 
%
