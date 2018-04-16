% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% y=Rmult(x)
%
% This function performs a matrix vector multiplication
%
%   y=R*x
%
% where R is the resolution matrix, R=G#*G.
%
% INPUT
%   x - the vector to multiply
%
% OUTPUT
%   y - the result of the multiplication
%
% GLOBAL
%   L - a roughening matrix
%
function y = Rmult(x);
global L;
global G;

% Some size parameters.
[r, s] = size(L);

% First, multiply G*x.
Gx = G * x;

% Now, solve a least squares problem to effectively multiply G#*G*x.
Gxpad = sparse([Gx; zeros(r, 1)]);
y = lsqr(@Gregmult, Gxpad, 1.0e-4, 5000);
