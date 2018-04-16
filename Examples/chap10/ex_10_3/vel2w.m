% Midpoint raypath calcuation function
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [vmid,wv]=vel2w(xmid,zmid,ip,kp,xn,zn,v);
%
% Generates a vector of raypath mid points velocities and the necessary
% functions (wv(:,i)) to calculate the Jacobian
%
%
function [vmid,wv]=vel2w(xmid,zmid,ip,kp,xn,zn,v)

vmid=zeros(length(xmid),1);
ip1=ip+1;
kp1=kp+1;

for k=1:length(xmid)

xf=(xmid(k)-xn(ip(k)))/(xn(ip1(k))-xn(ip(k)));
zf=(zmid(k)-zn(kp(k)))/(zn(kp1(k))-zn(kp(k)));
xf1=1.0-xf;
zf1=1.0-zf;

wv(k,1)=xf1*zf1;
wv(k,2)=xf*zf1;
wv(k,3)=xf1*zf;
wv(k,4)=xf*zf;

%  calculate velocity
vmid(k)=wv(k,1)*v(ip(k),kp(k))+wv(k,2)*v(ip1(k),kp(k))+wv(k,3)*v(ip(k),kp1(k))+wv(k,4)*v(ip1(k),kp1(k));

end
