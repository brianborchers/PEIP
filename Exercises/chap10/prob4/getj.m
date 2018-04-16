% Parameter Estimation and Inverse Problems, 2nd edition, 2011 
% by R. Aster, B. Borchers, C. Thurber
%
%
% [J,ttcal,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit)
%
% J is the Jacobian function of travel time partial derivatives with respece
%to nodes that bound the specified raypath
%rpstore contains the ray path coordinates and ttcal contains the
%calculated travel times on return
%calcrays=0 reinitializes the ray paths; calcrays~=1 uses the starting raypaths stored in the 3-d array rpinit.
%
function [J,ttcal,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit)

NSEG=PSCALE*2;
J=zeros(PSCALE*PSCALE,PSCALE*PSCALE);
ivec=2:NSEG;

rpstore=zeros(PSCALE^2,NSEG+1,2);
ttcal=zeros(PSCALE,PSCALE);

%  loop over sources and receivers
rpnum=0;
for j=1:PSCALE
xs=sc(j,1);
zs=sc(j,2);
for k=1:PSCALE
%keep track of the ray path number being evaluated
rpnum=rpnum+1;
xr=rc(k,1);
zr=rc(k,2);
%  set up initial path
dx=(xr-xs)/NSEG;
dz=(zr-zs)/NSEG;
xp=[xs:dx:xr];

zp=zs+[0:NSEG]*dz;

%initial rays not used?
if calcrays==0
rp=[xp;zp]';
else
rp(:,:)=rpinit(rpnum,:,:);
end

%  travel time modeling convergence loop
for it=1:NIT

%  now perturb path
rpnew=rp;

x2a=0.5*(rp(ivec+1,1)+rp(ivec-1,1));
z2a=0.5*(rp(ivec+1,2)+rp(ivec-1,2));
dxa=rp(ivec+1,1)-rp(ivec-1,1);
dza=rp(ivec+1,2)-rp(ivec-1,2);
dna=dxa.*dxa+dza.*dza;
rdx=dxa./sqrt(dna);
rdz=dza./sqrt(dna);
[cxul,czul]=cellfunc(x2a,xn,z2a,zn);
vmid=vel2(x2a,z2a,cxul,czul,xn,zn,v);
[vx,vz]=vel2d(x2a,z2a,cxul,czul,xn,zn,v);
% pseudo-bending calculations
vrd=vx.*rdx+vz.*rdz;
rvx=vx-vrd.*rdx;
rvz=vz-vrd.*rdz;
rvs=sqrt(rvx.*rvx+rvz.*rvz);

%  begin path segment loop
for i=1:NSEG-1

xxk=x2a(i);
zzk=z2a(i);

%  which node is to upper left of segment midpoint
if (rvs(i)~=0)
rvx(i)=rvx(i)/rvs(i);
rvz(i)=rvz(i)/rvs(i);
rcur=vmid(i)/rvs(i);
rtemp=rcur-sqrt(rcur*rcur-0.25*dna(i));

% compute the new points and distance of perturbations
% using convergence enhancement
xxk=x2a(i)+XFAC*rvx(i)*rtemp;
zzk=z2a(i)+XFAC*rvz(i)*rtemp;
% end of if rvs ne 0 section
end

%  store new path point coordinate
rpnew(i+1,1)=xxk;
rpnew(i+1,2)=zzk;

%  end of path segment loop
end

%  re-initialize rp ray paths
[n,m]=size(rp);

%check for ray path convergence
if norm(reshape(rp-rpnew,n*m,1))/norm(reshape(rp,n*m,1)) < CONV
break
end
rp=rpnew;

% end of ray tracing convergence loop
end
rpstore(rpnum,:,:)=rp;

jvec=[2:NSEG+1];
xmida=0.5*(rp(jvec,1)+rp(jvec-1,1));
zmida=0.5*(rp(jvec,2)+rp(jvec-1,2));
dxa=rp(jvec,1)-rp(jvec-1,1);
dza=rp(jvec,2)-rp(jvec-1,2);
ra=sqrt(dxa.*dxa+dza.*dza);
[cxul,czul]=cellfunc(xmida,xn,zmida,zn);
[vmid2w,wv]=vel2w(xmida,zmida,cxul,czul,xn,zn,v);

% loop over this ray path to compute travel time and derivatives
%contributions from the four surrounding nodes
nobs=j+PSCALE*(k-1);

for i=1:NSEG

%  which node is to upper left of segment midpoint
%construct G matrix indices
node1=cxul(i)+PSCALE*(czul(i)-1);
node2=cxul(i)+PSCALE*(czul(i)-1)+1;
node3=cxul(i)+PSCALE*czul(i);
node4=cxul(i)+PSCALE*czul(i)+1;

% update this (nobs) row of the G matrix
J(nobs,node1)=(ra(i)*wv(i,1))'+J(nobs,node1);
J(nobs,node2)=(ra(i)*wv(i,2))'+J(nobs,node2);
J(nobs,node3)=(ra(i)*wv(i,3))'+J(nobs,node3);
J(nobs,node4)=(ra(i)*wv(i,4))'+J(nobs,node4);

end

%get travel times
ttcal(j,k)=sum(ra./vmid2w);

%  end of source and receiver loops
end
end

