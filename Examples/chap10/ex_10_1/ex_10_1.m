% Example 10.1
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
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

%.load the particular realizations of tstor, ttobs and vtrue used in the book
load tstorenoise.mat

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

% the coordinates of the fast and slow centers
xfast = xn(3*PSCALE/4);
zfast = zn(3*PSCALE/4);
xslow = xn(PSCALE/2);
zslow = zn(PSCALE/2);

%construct the true model
vtrue=zeros(PSCALE,PSCALE);
vbackground=2900;
for i=1:PSCALE
  for j=1:PSCALE
    vtrue(i,j)=vbackground*...
      (1+0.10*exp(-.00004*((xn(i)-xfast)^2+(zn(j)-zfast)^2)))*...
      (1-0.15*exp(-.00004*((xn(i)-xslow)^2+(zn(j)-zslow)^2)));
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
bookfonts
colormap(gray)
caxis([2600 3200])
colorbar('FontSize',18);

%plot the ray paths for the true model
hold on
%forward problem bent-ray travel time modeling and plotting function
tstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,vtrue,sc,rc);
hold off
disp('Displaying the true model and ray paths (fig. 1)')
print -deps2 c10fcrossraypaths.eps


%size of model space
n=PSCALE^2;

%construct roughening matrix
L=zeros(n,n);
k=1;
for j=1:PSCALE
    for i=1:PSCALE;
        mtmp=zeros(PSCALE,PSCALE);
        % if possible add the non central entries
        if i>1;      mtmp(i-1,j  )=1; end;
        if j>1;      mtmp(i  ,j-1)=1; end;
        if i<PSCALE; mtmp(i+1,j  )=1; end;
        if j<PSCALE; mtmp(i  ,j+1)=1; end;
        L(k,:)=reshape(mtmp,1,n);
        k=k+1;
    end
end

%set values for L diagonal elements, taking into account edges and corners
for i=1:n
    L(i,i)=-sum(L(i,:));
end

LL=L'*L;

% the true slowness model as a vector
mtrue=reshape(1./vtrue,PSCALE*PSCALE,1);

%piter corresponds to a particular choice of alpha
piter=0;
piteralphas=10.^(-2:0.25:1.75);
NPITER=length(piteralphas);
[nn,mm]=size(vi);
vstore=zeros(NPITER,nn,mm);
alphasq=zeros(NPITER,1);
misfit=zeros(NPITER,MAXITER);
mnorm=zeros(NPITER,MAXITER);
mrms=zeros(NPITER,MAXITER);

disp('Displaying the model at each iteration and alpha (fig. 2)')
%display headings for the number we will displaying every iteration
disp('alpha, iteration, rmsnew, sqrt(chi^2)')

%loop over values of the regularization parameter
for piteralpha=piteralphas
  piter=piter+1;

  %starting model in each loop is constant velocity/slowness
  v=vbackground*ones(PSCALE,PSCALE);
  %convert velocity to slowness
  slow=1./v;

  %iteratively updated slowness model vector
  m=reshape(slow,PSCALE*PSCALE,1);

  vstartvec=reshape(v,PSCALE*PSCALE,1);

  % reset some GN state 
  calcrays=0;
  rpinit=0;
  rms=zeros(MAXITER,1);

  % start the iterative model loop
  for iter=1:MAXITER
    ttcal=zeros(PSCALE,PSCALE);

    %get the Jacobian matrix, J, and forward modeled travel times 
    %for the working model, v
    if iter==1
      [J,ttcal,rpstore]=...
          getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit);
    else
      J=K;
      ttcal=tttry;
    end

    %save raypaths to start the next iteration
    calcrays=1;
    rpinit=rpstore;

    %set the regularization tradeoff parameter here
    if (iter == 1)
      rparam=norm(J)*piteralpha;
      alphasq(piter)=rparam;
    end

    %calculate the travel-time residual vector for this iteration
    r=ttcal-ttobs;
    rvec=reshape(r,PSCALE*PSCALE,1);
    rms(iter,1)=norm(rvec);

    %GN, explicit regularization solution
    J1= J'*J+rparam*LL;
    rhs=-(J'*rvec+rparam*LL*m);
    dm=J1\rhs;

    %model update
    m=m+dm;
    slow=reshape(m,PSCALE,PSCALE);
    v=1./slow;

    %run the forward model to document variance reduction at this step and save
    %the Jacobian for future steps
    [K,tttry,rpstore]=getj(PSCALE,NIT,CONV,XFAC,xn,zn,v,sc,rc,calcrays,rpinit);
    rtry=ttobs-tttry;
    dtry=reshape(rtry,PSCALE*PSCALE,1);
    rmstry=norm(dtry);

    % display the alpha, iteration, rmsnew and chi^2 values
    disp([num2str(sqrt(alphasq(piter))) ' ' num2str(iter) ' ' num2str(rmstry)...
        ' ' num2str(NOISE*PSCALE)])

    % store some of this for later analysis
    misfit(piter,iter)=rmstry;
    mnorm(piter,iter)=norm(L*m);
    mrms(piter,iter)= norm((mtrue-m));

    % interpolate the current velocity structure for plotting
    vi=interp2(xnm,znm,v,xnim,znim,'cubic');


    %figure 2 shows the model for each iteration and alpha as the code runs
    figure(2)
    clf
    caxis([2300 3200])
    imagesc(xni,zni,vi')
    colormap('gray')
    colorbar('FontSize', 18);
    title(['\alpha = ',num2str(sqrt(alphasq(piter))),' iteration = '...
        ,num2str(iter),' rms residual = ',num2str(rmstry)]);
    bookfonts
    drawnow;

    %  end of GN iteration
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

disp('Displaying misfits by iteration (fig. 3)')

%this figure shows the seminorm statistics for each iteration and alpha value
figure(4)
clf
semilogy(1:iter,mnorm)
xlabel('Iteration')
ylabel('Solution Seminorm ||Lm||_2')
bookfonts

disp('Displaying the model semi-norms by iteration (fig. 4)')

%L-curve figure
figure(5)
clf
loglog(misfit(1:piter,MAXITER),mnorm(1:piter,MAXITER),'ok-')
xlabel('Residual Norm ||G(m)-d||_2')
ylabel('Solution Seminorm ||Lm||_2')
axis([.001 .1 .00001 .001])
hold on
% plot the discrepancy principle line
loglog([NOISE*PSCALE,NOISE*PSCALE],[.00001,.001],'k--')
bookfonts
% label each alpha
for i=1:2:piter
  text(misfit(i,MAXITER),mnorm(i,MAXITER),['   ',num2str(sqrt(alphasq(i)))]);
end
% label the discrepancy principle line
text(NOISE*PSCALE,.000015,'   \delta=0.008')
disp('Displaying the L-curve (fig. 5)')
print -depsc2 c10fcrossL.eps

%Models compared to true model in norm difference
figure(6)
clf
semilogx(sqrt(alphasq),mrms(1:piter,MAXITER),'ok-')
xlabel('\alpha')
ylabel('||m_{true}-m||_2')
bookfonts
hold on
% mark the selected alpha
semilogx(sqrt(alphasq(8)),mrms(8,MAXITER),'ok','MarkerSize',14)
% label the selected and end alphas
for i=[1,8,16]
  H=text(sqrt(alphasq(i)),mrms(i,MAXITER)+1.2e-5,...
      ['  ',num2str(sqrt(alphasq(i)),'%5.2f')]);
  set(H,'FontSize',12)
end
ylim([0 3.1e-4])
hold off
disp('Displaying alpha versus the model misfit (fig. 6)')
print -deps2 c10fcrossmodmisfit.eps

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
disp('Displaying the recovered model for each alpha (fig. 7)')
print -deps2 c10fcrossmodelsuite.eps

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
H=colorbar('FontSize', 18);
hold on
tstor=plotraypaths(PSCALE,NIT,CONV,XFAC,xn,zn,vtmp,sc,rc);
hold off
disp('Displaying the discrepancy principle selected model (fig. 8)')
print -deps2 c10fcrosstradeoff.eps
