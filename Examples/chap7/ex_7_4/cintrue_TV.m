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
function cin=cintrue_TV(t)

nt=size(t,2);

off=zeros(3,1);
m=off;
d=off;

% the peak times
off(1)=130.;
off(2)=150.;
off(3)=190.;

% change from times to samples
offs(1)=49;
offs(2)=60;
offs(3)=74;

% the peak heights
m(1)=0.7;
m(2)=0.3;
m(3)=0.5;

%the peak widths
d(1)=4;
d(2)=6;
d(3)=7;


% compute the true concentration
cin=zeros(1,nt);

for j=1:3
n1=offs(j)-d(j);
n2=offs(j)+d(j);
for i=n1:n2
	cin(i)=m(j);
end
end


