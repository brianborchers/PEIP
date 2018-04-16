%ray tracing subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% vmid=vel2(xmid,zmid,ip,kp,xn,zn,v);
%
%2-d velocity interpolation function
% INPUT
%   xmid - the x midpoint of the ray segment
%   zmid - the z midpoint of the ray segment
%   ip   - the first index to the cell of v containing xmid, zmid
%   kp   - the second index to the cell of v containing xmid, zmid
%   xn   - the x cell divisions
%   zn   - the z cell divisions
%   v    - the velocity structure
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

  % calculate velocity as a weighted combination of the four cells nearest the 
  % midpoint
  vmid(k)=wv1*v(ip(k),kp(k))+...
          wv2*v(ip1(k),kp(k))+...
          wv3*v(ip(k),kp1(k))+...
          wv4*v(ip1(k),kp1(k));
end
