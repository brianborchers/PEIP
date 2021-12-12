% Example 7.4
% This routine will return the true source
% history for the Skaggs and Kabala paper [reference 147 or 148].
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% function cin=cintrue(t)
%
% t = time array at which true solution is sought
% cin = true concentration vector
%
function cin=cintrue(t)

nt=size(t,2);

off=zeros(3,1);
m=off;
d=off;

% the peak times
off(1)=130.;
off(2)=150.;
off(3)=190.;

% the peak heights
m(1)=1.;
m(2)=0.3;
m(3)=0.5;

%the peak variances
d(1)=5;
d(2)=10;
d(3)=7;


% compute the sum of the described normal curves
cin=zeros(1,nt);
for i=1:3
	cin=cin+m(i)*exp(-(t-off(i)).^2/(2*d(i)^2));
end
