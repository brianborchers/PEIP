% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% p=mytcdf(x,nu)
%
% Computes the probability that a t random variable with nu degrees of
% freedom is less than or equal to x.  The computation is done by
% transformation to the standard normal, using a formula from Thistead,
% Elements of Statistical Computing, p334.
%
function p = mytcdf(x, nu)
s = 0.368 * (8 * nu + 3) / (2 * sqrt(nu^2*log(1+x^2/nu)));
Z = (1 - (2 * sqrt(1-exp(-s^2))) / (8 * nu + 3)) * sqrt(nu*log(1+x^2/nu));
p = phi(Z);
