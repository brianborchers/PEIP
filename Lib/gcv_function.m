% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% function g = gcv_function(alpha,s2,beta,delta0,mn)
%
% Auxiliary routine for GCV calculations
%
% INPUT
%   alpha -
%   s2     - the square of the gamma from the gsvd
%            function; these are in descending order
%   beta   - the projected data to fit
%   delta0 - the unavoidable misfit
%   mn     - the number of rows of U minus the number of columns of U
%            used as part of estimating the trace of I-GG#
%
% OUTPUT
%   g - || Gm_(alpha,L) - d ||^2 / (Tr(I - GG#)^2
function g = gcv_function(alpha, gamma2, beta)

f = (alpha^2) ./ (gamma2 + alpha^2);
length(f);
length(beta);
if length(f) > length(beta)
    f = f(1:length(beta));
else
    if length(beta) > length(f)
        beta = beta(end-length(f)+1:end);
    end
end
g = (norm(f.*beta)^2) / (sum(f))^2;
