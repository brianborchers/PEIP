% Example 10.4 subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% ufinal=forward(m,n,deltat,deltax,D,u0)
%
% Integrate the heat equation forward in time using the Crank-Nicolson
% implicit Euler method.
%
% The system of equations for each step in time is of the form
%
%    A*u^(k+1)=B*u^(k)
%
% where A and B are tridiagonal matrices.  
%
% The input parameters are:
%
%    n         Number of points in the x grid
%    m         Number of time steps.
%    deltat    Length of each time step. 
%    deltax    Grid spacing (deltax=1/(n+1))
%
% Outputs are:
%    u         The solution at time T=m*deltat
%    A,B       The Crank-Nicolson sparse matrices.
%
function [u,A,B]=forward(m,n,deltat,deltax,D,u0)
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
%
% Do the time steps.
%
u=u0;
for i=1:m
  u=A\(B*u);
end

