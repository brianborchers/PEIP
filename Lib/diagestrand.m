% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% d=diagestrand(matvecmult,n,s)
%
% This function uses the algorithm of Bekas et al. to estimate the diagonal
% of a square matrix, where the matrix is not given explicitly, but instead
% we have a function for computing matrix vector productions.
%
% Inputs:
%   matvecmult            The handle for a function that does multiplications
%   n                     The dimension of the matrix.
%   s                     The number of random test vectors to use.
%
% Output:
%   d                     The estimated diagonal of the matrix
%
function d = diagestrand(matvecmult, n, s)

% the sum of (A*v).*v and v.*v
t = zeros(n, 1);
q = t;

% generate s random test vectors
for k = 1:s
    v = randn(n, 1);
    t = t + matvecmult(v) .* v;
    q = q + (v .* v);
end
d = t ./ q;
