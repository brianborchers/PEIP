% Example 9.2.
% Computes the weighted jacobian of rawfunc.m
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% J=jac(p)
%
%
function J=jac(p)
global x;
global sigma;

m=length(p);
n=length(x);
J=zeros(m,n);
for j=1:n
    J(:,j)=(...
        [exp(p(2)*x(j));
        p(1)*x(j)*exp(p(2)*x(j));
        x(j)*exp(p(4)*x(j));
        p(3)*x(j)*x(j)*exp(p(4)*x(j))...
        ])/sigma(j);
end
J=J';
