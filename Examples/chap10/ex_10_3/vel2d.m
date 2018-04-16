% Derivative interpolation function for refreacted ray tomography
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [dvx,dvz]=vel2d(x2,z2,ip,kp,xn,zn,v);
%
%2-d derivative interpolation function
%
%
function [dvx,dvz]=vel2d(x2,z2,ip,kp,xn,zn,v)

dvx=zeros(length(x2),1);
dvz=zeros(length(x2),1);

ip1=ip+1;
kp1=kp+1;

for k=1:length(x2)

xd=xn(ip1(k))-xn(ip(k));
zd=zn(kp1(k))-zn(kp(k));

xf=(x2(k)-xn(ip(k)))/(xn(ip1(k))-xn(ip(k)));
zf=(z2(k)-zn(kp(k)))/(zn(kp1(k))-zn(kp(k)));
xf1=1.0-xf;
zf1=1.0-zf;
dvx(k)=(zf1*(v(ip1(k),kp(k))-v(ip(k),kp(k)))+zf*(v(ip1(k),kp1(k))-v(ip(k),kp1(k))))/xd;
dvz(k)=(xf1*(v(ip(k),kp1(k))-v(ip(k),kp(k)))+xf*(v(ip1(k),kp1(k))-v(ip1(k),kp(k))))/zd;

end
