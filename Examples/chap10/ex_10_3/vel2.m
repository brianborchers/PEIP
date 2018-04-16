% velocity interpolation function for refracted ray tomography
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% vmid=vel2(xmid,zmid,ip,kp,xn,zn,v);
%
%2-d velocity interpolation function
%
function vmid=vel2(xmid,zmid,ip,kp,xn,zn,v)

vmid=zeros(length(xmid),1);
ip1=ip+1;
kp1=kp+1;

for k=1:length(xmid)

xf=(xmid(k)-xn(ip(k)))/(xn(ip1(k))-xn(ip(k)));
zf=(zmid(k)-zn(kp(k)))/(zn(kp1(k))-zn(kp(k)));
xf1=1.0-xf;
zf1=1.0-zf;

wv1=xf1*zf1;
wv2=xf*zf1;
wv3=xf1*zf;
wv4=xf*zf;

%  calculate velocity
vmid(k)=wv1*v(ip(k),kp(k))+wv2*v(ip1(k),kp(k))+wv3*v(ip(k),kp1(k))+wv4*v(ip1(k),kp1(k));
end
