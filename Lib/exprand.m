% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% r=exprand(mu,m,n)
%
% generate an m by n array of exponentially distributed random variable
% realizations with expected value mu.  If n is not specified its
% default value is the same as m.  If neither m nor n is specified,
% then default values are 1.
%
function r = exprand(mu, m, n)
%
% Handle default values of m and n.
%
if (nargin == 1)
    m = 1;
    n = 1;
end
if (nargin == 2)
    n = m;
end
%
% Generate uniform random values, and apply the exponential inverse CDF.
%
r = -mu * log(rand(m, n));
