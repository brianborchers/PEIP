% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% lp=logprior(m)
%
% For this problem, our prior is uniform on m1=[0 2], m2=[-1 0], m3=[0 2],
% m4=[-1 0]
%
function lp=logprior(m)
if (m(1)>=0) && (m(1)<=2) && (m(2)>=-0.9) && (m(2)<=0) && (m(3)>=0) && (m(3)<=2) && (m(4)>=-0.9) && (m(4)<=0)
  lp=0;
else
  lp=-Inf;

end
