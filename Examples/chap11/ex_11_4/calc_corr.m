% Example 11.4 correlation function calculator
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
function c = calc_corr(x,laglen)
%
%c = calc_corr(x,laglen)
%
%returns the first laglen elements of the circular (normalized) crosscorrelation of the column vector x
%
c=zeros(laglen,1);
x=x-mean(x);
c(1)=dot(x,x);
for i=2:laglen+1
    c(i)=dot(x,circshift(x,i-1));
end
c=c/c(1);
c=[flipud(c(2:end));c];
