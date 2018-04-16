%ray tracing subroudine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [vmid,wv]=vel2w(xmid,zmid,ip,kp,xn,zn,v);
%
% Generates a vector of raypath mid points velocities and the necessary
% functions (wv(:,i)) to calculate the Jacobian
%
% INPUT
%   xmid - the x midpoint of the ray segment
%   zmid - the z midpoint of the ray segment
%   ip   - the first index to the cell of v containing xmid, zmid
%   kp   - the second index to the cell of v containing xmid, zmid
%   xn   - the x cell divisions
%   zn   - the z cell divisions
%   v    - the velocity structure
%
% OUTPUT
%   vmid - the interpolated velocity for the points
%   wv   - the factors needed to compute the Jacobian
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

  % store the velocities dependence on model elements
  wv(k,1)=xf1*zf1;
  wv(k,2)=xf*zf1;
  wv(k,3)=xf1*zf;
  wv(k,4)=xf*zf;

  %  calculate velocity
  vmid(k)=wv(k,1)*v(ip(k),kp(k))+...
          wv(k,2)*v(ip1(k),kp(k))+...
          wv(k,3)*v(ip(k),kp1(k))+...
          wv(k,4)*v(ip1(k),kp1(k));
end
