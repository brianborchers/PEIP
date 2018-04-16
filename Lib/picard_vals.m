% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [utd,utd_norm] = picard_vals(U,sm,d)
%
% return Picard plot parameters for subsequent plotting
%
% INPUT
%   U  - the U matrix from the SVD or GSVD
%   sm - singular values in decreasing order, or the lambdas divided by the mus
%        in decreasing order
%   d  - the data to fit
%
% OUTPUT
%   utd      - the columns of U transposed times d
%   utd_norm - utd./sm
%
function [utd, utd_norm] = picard_vals(U, sm, d)

k = length(sm);
for i = 1:k
    utd(i) = U(:, i)' * d;
    utd_norm(i) = utd(i) / sm(i);
end
