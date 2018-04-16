% travel time calculation routine for refracted ray tomography
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% ttstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc);
%
% ttstor is the matrix of travel times corresponding to the ray tracing
% for the source points sc, receiver points rc, and velocity model v
% indexed at xn and zn.  Size of the problem is PSCALE^2 x PSCALE^2
%
function ttstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc)
bookfonts

ttstor=zeros(PSCALE,PSCALE);

%loop over sources and receivers
for j=1:PSCALE
  xs=sc(j,1);
  zs=sc(j,2);
  for k=1:PSCALE
    xr=rc(k,1);
    zr=rc(k,2);

    nseg=PSCALE*2;
    %set up initial path
    dx=(xr-xs)/nseg;
    dz=(zr-zs)/nseg;
    xp=[xs:dx:xr];
    
    zp(1:nseg+1)=zs+[0:nseg]*dz;
    
    rp=[xp;zp]';
    
    %evaluate initial TT estimate (straight line path)
    tt=0;
    
    %loop over path segment midpoints
    for i=2:nseg+1
      xmid=0.5*(rp(i,1)+rp(i-1,1));
      zmid=0.5*(rp(i,2)+rp(i-1,2));  
      
      %  which node is to upper left of segment midpoint
      [cxul,czul]=cellfunc(xmid,xn,zmid,zn);
      vmid=vel2(xmid,zmid,cxul,czul,xn',zn',v);
      tt=tt+(sqrt(dx*dx+dz*dz))/vmid;
      
      %  end of path segment loop 
    end
    
    %  end of initial ray path section, save the result
    ttlast=tt;
    
    %  bending convergence loop
    for it=1:NIT
      
      tt=0;
      
      %  now perturb path
      rpnew=rp;
      
      %  loop over path segments
      for i=2:nseg
        x2=0.5*(rp(i+1,1)+rp(i-1,1));
        z2=0.5*(rp(i+1,2)+rp(i-1,2));
        xxk=x2;
        zzk=z2;
        
        %which node is to the upper left and above the segment midpoint
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
          xxk=x2+rvx*rtemp;
          zzk=z2+rvz*rtemp;
          
          % convergence enhancement
          xxk=XFAC*(xxk-x2)+x2;
          zzk=XFAC*(zzk-z2)+z2;
          
          % end of if rvs ne 0 section
        end
        
        %  store new path point coordinate
        rpnew(i,:)=[xxk,zzk];
        
        %  end of path segment loop
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
      
      %check for travel time convergence
      if abs(ttlast-tt)/ttlast < CONV
        break
      end
      
      %no convergence yet, save travel time and reiterate
      ttlast=tt;
      
    end
    
    %store tt
    ttstor(j,k)=tt;
    
    plot(rp(:,1),rp(:,2),'k--')
    axis('square')
    axis('equal')
    axis ij
    axis([min(xn),max(xn),min(xn),max(xn)])
    xlabel('m')
    ylabel('m')
    
    %  end of source and receiver loops
  end
end
