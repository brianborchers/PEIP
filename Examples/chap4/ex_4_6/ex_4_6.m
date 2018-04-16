% Example 4.6
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% (originally suggested by Gary Pavlis, Indiana University).

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% load the saved data
load vsp

% build the system model
% the number of data points
M=length(t);
% Number of model parameters relative to number of observations (k is an
% integer greater than one)
k=1;
N=k*M;
dd=1000/N;

%System matrix
G=zeros(M,N);
for i=1:M
  G(i,1:k*i) = ones(1,k*i)*dd;
end
%G=tril(ones(N,MM))*dd;

%first order Tikh regularization setup
L1 = get_l_rough(N,1);
[U1,V1,X1,Lam1,M1]=gsvd(G,L1);

% get the points ans solutions for the first order TGSVD L-curve
[rho1,eta1,reg_param1,m1s]=l_curve_tgsvd(U1,t,X1,Lam1,G,L1);

% store where the corner is (from visual inspection)
ireg_corner1=8;
rho_corner1=rho1(ireg_corner1);
eta_corner1=eta1(ireg_corner1);
disp(['1st order reg corner is:  ',num2str(ireg_corner1)]);


% plot the L-curve for TGSVD regularization
% note that the k index starts at zero
figure(1)
clf
loglog(rho1,eta1,'ko', rho1,eta1,'k-.');
xlim([1e-4 1e-2])
ylim([6e-6 2e-4])
xlabel('Residual Norm ||Gm-d||_2')
ylabel('Solution Seminorm ||Lm||_2')
bookfonts

for i=[10 20 40],
  H=text(rho1(i),eta1(i),['    ',num2str(i)]);
  set(H,'FontSize',18);
end
hold on
%plot the corner on the L-curve
H=loglog(rho_corner1,eta_corner1,'ko');
set(H,'markersize',16)
print -deps2 c4flcurve1tsvd.eps

disp('Displaying TGSVD first order L-curve (fig. 1)');

%second order Tikh regularization setup
L2 = get_l_rough(N,2);
[U2,V2,X2,Lam2,M2]=gsvd(G,L2);

% get the points and solutions for the second order TGSVD L-curve
[rho2,eta2,reg_param,m2s]=l_curve_tgsvd(U2,t,X2,Lam2,G,L2);

% store where the corner is (from visual inspection)
% note that the k index starts at zero
ireg_corner2=7;
rho_corner2=rho2(ireg_corner2);
eta_corner2=eta2(ireg_corner2);

disp(['2nd order reg corner is:  ',num2str(ireg_corner2)]);

%plot the second order L-curve based on the TGSVD method
figure(2)
clf
loglog(rho2,eta2,'ko',rho2,eta2,'k-.');
xlim([1e-4 1e-2]);
ylim([1e-7 1e-3]);
xlabel('Residual Norm ||Gm - d||_2');
ylabel('Solution Seminorm ||Lm||_2');
bookfonts

for i=[10 20 40],
  HH=text(rho2(i),eta2(i),['    ',num2str(i)]);
  set(HH,'fontsize',18);
end
hold on
%plot the corner on the L-curve
H=loglog(rho_corner2,eta_corner2,'ko');
set(H,'markersize',16)
print -deps2 c4flcurve2tsvd.eps

disp('Displaying TGSVD second order L-curve (fig. 2)');

% get the models
m1=m1s(:,ireg_corner1);
m2=m2s(:,ireg_corner2);

% plot the L1 model results
figure(3)
clf
plotconst(m1*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'--')
xlabel('Depth (m)')
ylabel('Slowness (s/km)')
bookfonts
ylim([0.24 0.34]);
print -deps2 c4fmtik1tsvd.eps

disp('Displaying first order Tikhonov model (fig. 3)');

% plot the L2 model results
figure(4)
clf
plotconst(m2*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'k--')
xlabel('Depth (m)')
ylabel('Slowness (s/km)')
bookfonts
ylim([0.24 0.34]);
print -deps2 c4fmtik2tsvd.eps

disp('Displaying second order Tikhonov model (fig. 4)');

%compute and display the model misfits relative to a highly-sampled smooth
%model evaluated at the midpoints of each constant slowness model segment
dinterp=dd-dd/2:dd:maxdepth-dd/2;
strueresamp=slowfunctionsmooth(dinterp)';

disp(['1st order reguarization, model 2-norm misfit from true model (s/km): ',num2str(1000*norm(m1-strueresamp))]);
disp(['2nd order regularization, model 2-norm misfit from true model (s/km): ',num2str(1000*norm(m2-strueresamp))]);
