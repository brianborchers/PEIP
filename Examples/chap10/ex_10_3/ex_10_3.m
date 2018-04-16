%Example 10.3 
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
clear

%number of iterations in ray bending
NIT=50;

%scale of problem (m)
SCALE=1600;

%scale of problem (receivers, sources and square root of number of velocity nodes)
PSCALE=8;

% convergence enhancement factor for subsequent raybending subroutines
XFAC=1.2;

%Number of GN linearization steps
MAXITER=5;

%Travel-time convergence parameter
CONV=1e-4;

%noise level (s)
NOISE = 0.001;

%  source and receiver coordinate matrices
scz=linspace(0,SCALE,PSCALE);
scx=zeros(size(scz));
sc=[scx;scz]';
rcz=linspace(0,SCALE,PSCALE);
rcx=SCALE*ones(size(scz));
rc=[rcx;rcz]';

%  velocity model node positions (xn,zn)
xn=linspace(-150,SCALE+150,PSCALE)';
zn=xn;

%construct the true model
%x,y indices
vtrue=zeros(PSCALE,PSCALE);
vbackground=2900;

%checkerboard vtrue model for resolution test
vtrue=vbackground*ones(size(vtrue));
for i=2:PSCALE-1
    for j=2:PSCALE-1
        vtrue(i,j)=vbackground*(1+(-1)^(i*j)*0.10);
    end
end

%linerly interpolate the velocity field for smooth plotting
[xnm,znm]=meshgrid(xn,zn');
xni=(-150:10:SCALE+150)';
zni=(-150:10:SCALE+150)';
[xnim,znim]=meshgrid(xni,zni);
vi=interp2(xnm,znm,vtrue,xnim,znim,'cubic');

%plot velocity structure and ray paths
figure(1)
clf
%plot the transpose to show (row, column) view of model
imagesc(xni,zni,vi')
colormap(gray)
bookfonts
caxis([2600 3200])
H=colorbar;
set(H,'FontSize',18);

%plot the ray paths for the true model
hold on

%forward problem bent-ray travel time modeling and plotting function
tstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,vtrue,sc,rc);
hold off

%print out resolution figure
print -deps2 c10fchecknlres1.eps

%size of model space
n=PSCALE^2;

%construct roughening matrix
L=zeros(n,n);
k=1;
for j=1:PSCALE
    for i=1:PSCALE;
        mtmp=zeros(PSCALE,PSCALE);
        if i>1; mtmp(i-1,j)=1; end;
        if j>1; mtmp(i,j-1)=1; end;
        if i<PSCALE; mtmp(i+1,j)=1; end;
        if j<PSCALE; mtmp(i,j+1)=1; end;
        L(k,:)=reshape(mtmp,n,1);
        k=k+1;
    end
end

%set values for L diagonal elements, taking into account edges and corners
for i=1:n
    L(i,i)=-sum(L(i,:));
end

%construct a first--order regularization matrix, if desired.
% L1=zeros(n,n);
% k=1;
% for j=1:PSCALE
%     for i=1:PSCALE;
%         mtmp=zeros(PSCALE,PSCALE);
%         if i>1; mtmp(i-1,j)=1; end;
%         if j>1; mtmp(i,j-1)=1; end;
%         L1(k,:)=reshape(mtmp,n,1);
%         k=k+1;
%     end
% end
% for i=1:n
%     L1(i,i)=-sum(L1(i,:));
% end


% constuct a hybrid regularization matrix, if desired
% L=L+L1;

%piter corresponds to a particular choice of alpha
piter=0;
piterexprange=(-2:0.25:1.75);
NPITER=length(piterexprange);
[nn,mm]=size(vi);
vstore=zeros(NPITER,nn,mm);
alphasq=zeros(NPITER,1);
misfit=zeros(NPITER,MAXITER);
mnorm=zeros(NPITER,MAXITER);
mrms=zeros(NPITER,MAXITER);

%loop over values of the regularization parameter
for piterexp=piterexprange
    
piter=piter+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now evaluate the travel time data for this true model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%starting model in each loop is constant velocity/slowness
v=vbackground*ones(PSCALE,PSCALE);
%convert velocity to slowness
slow=1./v;

%iteratively updated slowness model vector
m=reshape(slow,PSCALE*PSCALE,1);

wv=zeros(4);

%noise-free data and true velocity model
%load tstor.mat

%noist-free data for checkerboard resoution Test
ttobs=tstor;

%generate another realization of observed travel times with noise
%ttobs=tstor+NOISE*randn(size(tstor));
%save tstorenoise.mat ttobs

%...or load a particular realization used in the textbook
%load tstorenoise.mat

mtrue=reshape(1./vtrue,PSCALE*PSCALE,1);

vstartvec=reshape(v,PSCALE*PSCALE,1);

%path segments
NSEG=PSCALE*2;

% start of iterative model loop
calcrays=0;
rpinit=0;

LL=L'*L;

for iter=1:MAXITER

%ttcal=zeros(PSCALE,PSCALE);

%get the Jacobian matrix, J, and forward modeled travel times for the working model, v
 if iter==1
     [J,ttcal,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit);
 else
 J=K;
 ttcal=tttry;
 end

%save raypaths to start the next iteration
  calcrays=1;
  rpinit=rpstore;

%set the regularization tradeoff parameter here

  if (iter == 1)
    rparam=norm(J)*10^piterexp;
  end

  alphasq(piter)=rparam;

%calculate the travel-time residual vector for this iteration
  rms=zeros(MAXITER,1);
  r=ttcal-ttobs;
  rvec=reshape(r,PSCALE*PSCALE,1);
  rms(iter,1)=sqrt(rvec'*rvec);

%GN, explicit regularization solution
  J1= J'*J+rparam*LL;
  rhs=-(J'*rvec+rparam*LL*m);
  dm=J1\rhs;

%model update
  m=m+dm;
  slow=reshape(m,PSCALE,PSCALE);
  mnslow=mean(mean(slow));
  v=1./slow;

%run the forward model to document variance reduction at this step and save
%the Jacobian for future steps if needed
  [K,tttry,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit);
  rtry=ttobs-tttry;
  dtry=reshape(rtry,PSCALE*PSCALE,1);
  rmstry=sqrt(dtry'*dtry);
  disp('alpha, iteration, rmsnew, sqrt(chi^2)')
  [sqrt(alphasq(piter)) iter rmstry NOISE*PSCALE]
  misfit(piter,iter)=rmstry;
  mnorm(piter,iter)=norm(L*m);

%rms difference wrt true model
  mrms(piter,iter)= norm((mtrue-m));

%figure 2 shows the separate model for each iteration as the code runs, for different regularization parameters (alpha).
  figure(2)
  clf
%interpolate the velocity field for plotting
  vi=interp2(xnm,znm,v,xnim,znim,'cubic');
  caxis([2300 3200])
  imagesc(xni,zni,vi')
  bookfonts
  colormap('gray')
  H=colorbar;
  set(H,'FontSize',18);
  title(['\alpha = ',num2str(sqrt(alphasq(piter))),' iteration = ',num2str(iter),' rms residual = ',num2str(rmstry)]);
  drawnow;
  pause(1);

%  end of tomography inversion loop
end

%store this model
vstore(piter,:,:)=vi;

% end of regularization parameter loop
end

%show results for the various regularization parameters

%this figure shows the residual stats for each iteration and alpha value
figure(3)
clf
plot(1:iter,misfit)
xlabel('Iteration')
ylabel('Residual Norm ||G(m)-d||_2')
bookfonts

%this figure shows the seminorm statistics for each iteration and alpha value
figure(4)
clf
semilogy(1:iter,mnorm)
xlabel('Iteration')
ylabel('Solution Seminorm, ||Lm||_2')
bookfonts

%L-curve figure
figure(5)
clf
loglog(misfit(1:piter,MAXITER),mnorm(1:piter,MAXITER),'ok-')
xlabel('Residual Norm ||G(m)-d||_2')
ylabel('Solution Seminorm ||Lm||_2')
axis([.001 .1 .00001 .001])
hold on
loglog([NOISE*PSCALE,NOISE*PSCALE],[.00001,.001],'k--')
bookfonts
for i=1:2:piter
text(misfit(i,MAXITER),mnorm(i,MAXITER),['   ',num2str(sqrt(alphasq(i)))]);
end
text(NOISE*PSCALE,.000015,'   \delta=0.008')
hold off

%Models compared to true model in norm difference
figure(6)
clf
semilogx(sqrt(alphasq),mrms(1:piter,MAXITER),'ok-')
xlabel('\alpha')
ylabel('||m_{true}-m||_2')
hold on
semilogx(sqrt(alphasq(8)),mrms(8,MAXITER),'ok','MarkerSize',14)
bookfonts
for i=[1,8,16]
H=text(sqrt(alphasq(i)),mrms(i,MAXITER)+1.2e-5,['  ',num2str(sqrt(alphasq(i)),'%5.2f')]);
set(H,'FontSize',12)
end
ylim([0 3.1e-4])
hold off

%Suite of models figure (color on cover of book)
figure(7)
clf
for i=1:piter
subplot(4,4,i)
vtmp1=vstore(i,:,:);
vtmp(:,:)=vtmp1;
imagesc(vtmp')
caxis([2600 3200])
axis square
set(gca,'xticklabel','')
set(gca,'yticklabel','')
title(['\alpha = ',num2str(sqrt(alphasq(i)))])
colormap('gray')
end

%print out resolution figure
print -deps2 c10fchecknlres2.eps

%"best" model figure from discrepancy principle and ray paths
figure(8)
clf
vtmp1=vstore(8,:,:);
vtmp(:,:)=vtmp1;
imagesc(xni,zni,vtmp')
bookfonts
caxis([2600 3200])
axis square
xlabel('m')
ylabel('m')
colormap('gray')
H=colorbar;
set(H,'FontSize',18);
hold on
tstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,vtmp,sc,rc);
hold off

