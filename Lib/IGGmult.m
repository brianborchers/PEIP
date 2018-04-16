% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% This function performs the matrix vector multiplication
%
%   (I-G*G#)*x
%
% required by diagestrand() for use in computing the GCV function.
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
%
function y = IGGmult(x)
global L G;

% Some size parameters.
[r, ~] = size(L);

% First, multiply G#*x by solving the least squares problem with LSQR.
% Note that we don't require a very high accuracy solution here.
xpad = [x; zeros(r, 1)];
y = lsqr(@Gregmult, xpad, 1.0e-4, 1500);

% Now, multiply G*G#*x.
y = G * y;

% Finally, let y=x-y to get y=(I-G*G#)*x;
y = x - y;
