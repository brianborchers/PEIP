% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% Function r=randt(nu,m,n) generates an m by n matrix of random numbers
% with a t distribution with nu degrees of freedom.
% (Uses function tinv.m.)
%
% Input Parameters:
%   nu - degrees of freedom of t distribution
%   m - row dimension of returned matrix
%   n - column dimension of returned matrix
%
% Output Parameters:
%   r - m x n matrix of random t distribution values.
%
function r = randt(nu, m, n)

% If m,n not specified, then use 1,1.
if (nargin == 1),
    m = 1;
    n = 1;
end

% Initialize matrices
U = rand(m, n);
r = zeros(m, n);

% Populate matrix and return
for i = 1:m
    for j = 1:n
        r(i, j) = tinv(U(i, j), nu);
    end
end
