% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% [v0,num,denom,normgmd,normlm]=gcv(a)
%
% Computes the generalized cross validation function for a given value, a, of
% the regularization parameter, using the stochastic approach of Bekas et al.
% to estimate the required trace.
%
% INPUT
%   a - the value of the regularization parameter to use
%
% OUTPUT
%   v0      - the value of the GCV function
%   num     - the numerator of v
%   denom   - the denominator of v
%   normgmd - the norm of the misfit norm(G*m-d)
%   normlm  - the model seminorm norm(L*m)
%
% GLOBAL
%   G     - the system matrix
%   L     - the regularization matrix
%   alpha - the regularization parameter (set to a on entry)
%   d     - the goal data
%
function [v0, num, denom, normgmd, normlm] = gcv(a)

% Global variables used within this function.
global alpha;
global G;
global L;
global d;

% Set alpha.
alpha = a;

% Get the size of G and L.
[m, n] = size(G);
[r, s] = size(L);

% First, compute the regularized solution.
dpad = sparse([d; zeros(r, 1)]);
mreg = lsqr(@Gregmult, dpad, [], 5000);

% Compute norm(G*mreg-d) and norm(L*mreg) so they can be returned
normgmd = norm(G*mreg-d);
normlm = norm(L*mreg);

% The numerator is norm(G*mreg-d)^2;
num = norm(G*mreg-d)^2;

% Now, work on the denominator.
numvecs = 256;
diag = diagestrand(@IGGmult, m, numvecs);

% Now, compute the denominator.
denom = sum(diag)^2;

% Put the numerator and denominator together to get v0.
v0 = m * num / denom;
