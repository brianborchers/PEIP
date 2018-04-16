%ray tracing subroutine
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% ttstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc);
%
% INPUT
%   PSCALE - the length of a side of v
%   NIT    - the maximum number of iterations to perform determing the path
%   CONV   - the relative change in travel time that is acceptable for 
%            convergence
%   XFAC   - a factor used to control convergance
%   xn     - the x positions of the nodes in v
%   zn     - the z positions of the nodes in v
%   v      - the seismic velocity grid (PSCALE by PSCALE matrix)
%   sc     - the coordinates of the seismic sources
%   rc     - the coordinates of the seismic receivers
%
% OUTPUT
%   ttstor - the travel times between each source and each receive
%
% ttstor is the matrix of travel times corresponding to the ray tracing
% for the source points sc, receiver points rc, and velocity model v
% indexed at xn and zn.  Size of the problem is PSCALE x PSCALE
%
% This function also adds the ray paths to the current plot.
%
function ttstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc)

%make sure the current plot has the correct fonts
bookfonts

ttstor=zeros(PSCALE,PSCALE);

% the number of segments a ray path will have
nseg=PSCALE*2;

% loop over sources
for j=1:PSCALE
  xs=sc(j,1);
  zs=sc(j,2);
  % loop over receivers
  for k=1:PSCALE
    xr=rc(k,1);
    zr=rc(k,2);

    % get the evenly split difference in each direction
    dx=(xr-xs)/nseg;
    dz=(zr-zs)/nseg;

    % set up initial path as a strait line
    xp=[xs:dx:xr];
    zp=zs+[0:nseg]*dz;
    
    % the initial path matrix
    rp=[xp;zp]';
    
    % evaluate initial travel time (TT) estimate
    tt=0;
    
    % estimate the travel time by summing the time for each segment
    for i=1:nseg
      xmid=0.5*(rp(i+1,1)+rp(i,1));
      zmid=0.5*(rp(i+1,2)+rp(i,2));  
      
      % which node is to upper left of segment midpoint
      [cxul,czul]=cellfunc(xmid,xn,zmid,zn);
      vmid=vel2(xmid,zmid,cxul,czul,xn',zn',v);
      tt=tt+(sqrt(dx*dx+dz*dz))/vmid;
    end
    
    % end of initial travel time section, save the result
    ttlast=tt;
    
    %  bending convergence loop
    for it=1:NIT
      tt=0;
      
      % compute a new path by modifying this one
      rpnew=rp;
      
      % loop over internal path points
      for i=2:nseg
        %get the midpoint between the previous and next points
        x2=0.5*(rp(i+1,1)+rp(i-1,1));
        z2=0.5*(rp(i+1,2)+rp(i-1,2));
        xxk=x2;
        zzk=z2;
        
        %which node is to the upper left of the strait line point
        [cxul, czul]=cellfunc(x2,xn,z2,zn);
        vmid=vel2(x2,z2,cxul,czul,xn,zn,v);
        [vx,vz]=vel2d(x2,z2,cxul,czul,xn,zn,v);
        
        % pseudo-bending calculations
        dx=rp(i+1,1)-rp(i-1,1);
        dz=rp(i+1,2)-rp(i-1,2);
        dn=dx*dx+dz*dz;
        ddn=sqrt(dn);
        rdx=dx/ddn;
        rdz=dz/ddn;
        vrd=vx*rdx+vz*rdz;
        rvx=vx-vrd*rdx;
        rvz=vz-vrd*rdz;
        rvs=sqrt(rvx*rvx+rvz*rvz);
        
        if (rvs~=0)
          rvx=rvx/rvs;
          rvz=rvz/rvs;
          rcur=vmid/rvs;
          rtemp=rcur-sqrt(rcur*rcur-0.25*dn);
          
          % compute the new points and distance of perturbations
          % using convergence enhancement
          xxk=x2+XFAC*rvx*rtemp;
          zzk=z2+XFAC*rvz*rtemp;
        end
        
        % store path point at the newly computed coordinate
        rpnew(i,:)=[xxk,zzk];
        
        % end of path interior point loop
      end
      
      %  re-initialize rp
      rp=rpnew;
      
      % loop over new path to compute travel time
      
      for i=2:nseg+1
        xmid=0.5*(rp(i,1)+rp(i-1,1));
        zmid=0.5*(rp(i,2)+rp(i-1,2));
        
        %which node is to the upper left and above the segment midpoint
        [cxul, czul]=cellfunc(xmid,xn,zmid,zn);
        vmid=vel2(xmid,zmid,cxul,czul,xn,zn,v);
        dx=rp(i,1)-rp(i-1,1);
        dz=rp(i,2)-rp(i-1,2);
        tt=tt+(sqrt(dx*dx+dz*dz))/vmid;
      end
      
      %check this ray's travel time for convergence
      if abs(ttlast-tt)/ttlast < CONV
        break
      end
      
      %no convergence yet, save travel time and reiterate
      ttlast=tt;
    end
    
    %store tt
    ttstor(j,k)=tt;
    
    % add this ray to the plot
    plot(rp(:,1),rp(:,2),'k--')
    axis('square', 'equal', 'ij')
    axis([min(xn),max(xn),min(xn),max(xn)])
    xlabel('m')
    ylabel('m')
  end
end
