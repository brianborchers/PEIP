% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%  y=Gregmult(x, trans)
%
% This computes the products
%
%   y=[G; alpha*L]*x
%
% and
%
%  y=[G' alpha*L')*x
%
% Required by LSQR to solve the least squares problem
%
% min norm(G*m-d)^2+norm(L*m)^2
%
% INPUT
%   trans - either 'notransp' or anything else
%   x     - a vector either with the same number of rows G and L have columns
%           or has the same number of rows as G plus L if trans is not notransp
%
% OUTPUT
%   y - the result of the multiplication as a column vector
%
% GLOBAL
%   G     - the system matrix
%   L     - the regularization matrix
%   GT    - the transpose of the system matrix
%   LT    - the transpose of the regularization matrix
%   alpha - the regularization parameter
%
%
function y = Gregmult(x, trans)

% Global variables hold G, L, and their transposes, plus alpha.
global G;
global GT;
global L;
global LT;
global alpha;

% See if LSQR wants to do a regular multiplication or a transposed
% multiplication.
%
if (strcmp(trans, 'notransp'))
    % A regular multiplication.
    y = [G * x; alpha * (L * x)];
else
    % A transposed multiplication.
    [p, q] = size(G);
    [r, s] = size(L);
    y = GT * x(1:p) + alpha * (LT * x(p+1:p+r));
end
