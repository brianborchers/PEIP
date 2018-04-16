%ray tracing subrouting
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [J,ttcal,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit)
%
% INPUT
%   PSCALE   - the length of a side of v
%   NIT      - the maximum number of iterations to perform determing the path
%   CONV     - the relative change in travel time that is acceptable for 
%              convergence
%   XFAC     - a factor used to control convergance
%   xn       - the x positions of the nodes in v
%   zn       - the z positions of the nodes in v
%   v        - the seismic velocity grid (PSCALE by PSCALE matrix)
%   sc       - the coordinates of the seismic sources
%   rc       - the coordinates of the seismic receivers
%   calcrays - if 0 reinitialize the ray paths otherwise start using the ray 
%              paths in the 3-d array rpinit.
% OUTPUT
%   J       - the Jacobian function of travel time partial derivatives with 
%             respect to nodes that bound the specified raypath
%   ttcal   - contains the calculated travel times on return
%   rpstore - contains the ray path coordinates and 
%
function [J,ttcal,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit)

% the number of segments a ray path will have
NSEG=PSCALE*2;

% an empty array to hold the jacobian
J=zeros(PSCALE^2,PSCALE^2);

% a vector of the internal points on a ray path
ivec=2:NSEG;

% an array to hold the points on each ray path
rpstore=zeros(PSCALE^2,NSEG+1,2);
% an array that holds the travel times calculated
ttcal=zeros(PSCALE,PSCALE);

% loop over source
rpnum=0;
for j=1:PSCALE
  xs=sc(j,1);
  zs=sc(j,2);
  % loop over receivers
  for k=1:PSCALE
    xr=rc(k,1);
    zr=rc(k,2);

    %keep track of the ray path number being evaluated
    rpnum=rpnum+1;

    % get the evenly split difference in each direction
    dx=(xr-xs)/NSEG;
    dz=(zr-zs)/NSEG;

    % set up initial path as a strait line
    xp=[xs:dx:xr];
    zp=zs+[0:NSEG]*dz;

    % start with the strait line ray path if we were told to
    if calcrays==0
      rp=[xp;zp]';
    else
      rp(:,:)=rpinit(rpnum,:,:);
    end

    % ray path convergence loop
    for it=1:NIT
      % compute a new path by modifying this one
      rpnew=rp;

      % let every interior point be the midpoint of its neighbors
      x2a=0.5*(rp(ivec+1,1)+rp(ivec-1,1));
      z2a=0.5*(rp(ivec+1,2)+rp(ivec-1,2));

      % compute all of the deltas 
      dxa=rp(ivec+1,1)-rp(ivec-1,1);
      dza=rp(ivec+1,2)-rp(ivec-1,2);

      dna=dxa.*dxa+dza.*dza;
      rdx=dxa./sqrt(dna);
      rdz=dza./sqrt(dna);

      % determine velocity structure information for each point
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
        end

        %  store new path point coordinate
        rpnew(i+1,1)=xxk;
        rpnew(i+1,2)=zzk;

        %  end of path segment loop
      end

      [n,m]=size(rp);
      %check for ray path convergence
      if norm(reshape(rp-rpnew,n*m,1))/norm(reshape(rp,n*m,1)) < CONV
        rp=rpnew;
        break
      end
      rp=rpnew;

      % end of ray tracing convergence loop
    end
    rpstore(rpnum,:,:)=rp;

    % for every line segment
    jvec=1:NSEG;
    % get the midpoints
    xmida=0.5*(rp(jvec+1,1)+rp(jvec,1));
    zmida=0.5*(rp(jvec+1,2)+rp(jvec,2));
    % get the lengths
    dxa=rp(jvec+1,1)-rp(jvec,1);
    dza=rp(jvec+1,2)-rp(jvec,2);
    ra=sqrt(dxa.*dxa+dza.*dza);
    % get the reference cell
    [cxul,czul]=cellfunc(xmida,xn,zmida,zn);
    % get the midpoint velocities and model parameter dependence information
    [vmid2w,wv]=vel2w(xmida,zmida,cxul,czul,xn,zn,v);

    % the row of the Jacobian to update
    nobs=j+PSCALE*(k-1);

    % loop over this ray path to compute travel time and derivatives
    % contributions from the four surrounding nodes
    for i=1:NSEG
      % get the indices of the four influencing model nodes
      node1=cxul(i)+PSCALE*(czul(i)-1);
      node2=cxul(i)+PSCALE*(czul(i)-1)+1;
      node3=cxul(i)+PSCALE*czul(i);
      node4=cxul(i)+PSCALE*czul(i)+1;

      % update this row of the J matrix for each influencing model node
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
