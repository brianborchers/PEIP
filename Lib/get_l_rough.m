% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% L = get_l_rough(n,deg)
%
% returns a 1D differentiating matrix operating on a series with n points.
%
% INPUT
%   n   - the number of model points
%   deg - the order of the derivative to approximate
%
% OUTPUT
%   L - the discrete differentiation matrix
%
function L = get_l_rough(n, deg)
if deg < 0 | floor(deg) ~= deg
    disp('degree must be a non-negative integer');
    return
end

if deg == 0
    L = eye(n);
    return
end

% let df approximate the first derivative
df = [-1, 1, zeros(1, deg-1)];
for i = 2:deg
    % take the difference of the lower order derivative and itself shifted left
    % to get a derivative one order higher
    df = [0, df(1:deg)] - [df(1:deg), 0];
end

dn = n - deg;
L = sparse(n-deg, n);

% add the ith element of df to L i-1 elements to the right of the diagonal
for i = 1:deg + 1
    L = L + sparse(1:n-deg, [1:dn]+i-1, df(i)*ones(1, dn), dn, n);
end
L = full(L);
