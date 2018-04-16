% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% acceptance probability function
%
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
function z=getlnprobacc(m1,c,x,y,sigma)
%calculate log likelihoods
ll1=(-1/2)*sum((y-fun(c,x)).^2./sigma);
ll2=(-1/2)*sum((y-fun(m1,x)).^2./sigma);
z=min(0,ll1-ll2);
return;
