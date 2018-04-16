% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% L=lcoeff(n)
%
% Returns an (n+1) by (n+1) matrix.  Column i gives the
% coefficients of the i-1 th Legendre polynomial in a form suitable for
% use with polyval.
%
function L = lcoeff(n)
if (n < 3)
    L = [];
    return
end
L = zeros(n+3, n+1);
L(n+1, 1) = 1;
L(n, 2) = 1;
for j = 3:n + 1
    for i = 1:n + 1
        L(i, j) = L(i, j) + (2 * j - 3) * L(i+1, j-1) / (j - 1) - (j - 2) * L(i, j-2) / (j - 1);
    end
end
L = L(1:(n + 1), 1:(n + 1));
