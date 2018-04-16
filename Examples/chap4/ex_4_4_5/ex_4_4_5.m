% Examples 4.4 and 4.5
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

load vsp.mat

k=1;
%Number of model parameters (set k equal to an integer >=1 to explore the
%n>m case in this example)
N=k*M;

dd=20/k;

% generate the true smooth model
strue=slowfunctionsmooth(depth);

% create the noisy observed data
%
% dobs is a vector containing the depths of the different observations.
% dd is difference in depth between successive observations.
% t contains the noisy observations.
%
%forward model borehole arrival time data using function 'slowfunctionsmooth'
i=0;

%System matrix
G=zeros(M,N);
for i=1:M
  G(i,1:k*i) = ones(1,k*i)*dd;
end

% Example 4.4
figure(1)
clf
plot((0:length(strue)-1),strue*1000,'k');
xlabel('Depth (m)');
ylabel('True Slowness (s/km)');
bookfonts
print -deps2 c4fvspmod.eps

disp('Displaying the true slowness model (fig. 1)')

% Compute the least squares solution of G
mls=G\t;

%plot the least-squares solution and confidence intervals
figure(2)
clf
plotconst(mls*1000,0,maxdepth);
hold on
plotconstc(mls*1000-1.96*1000*0.0002*sqrt(diag(inv(G'*G))),0,maxdepth,'k--');
plotconstc(mls*1000+1.96*1000*0.0002*sqrt(diag(inv(G'*G))),0,maxdepth,'k--');
xlabel('Depth (m)');
ylabel('Slowness (s/km)');
bookfonts
print -deps2 c4fml2.eps

disp('Displaying the least squares model (fig. 2)')

%apply first-order Tikhonov regularization
L1 = get_l_rough(N,1);
[U1,V1,X1,LAM1,MU1]=gsvd(G,L1);

% apply the L curve criteria to the first-order regularization problem
[rho1,eta1,reg_param1,m1s]=l_curve_tikh_gsvd(U1,t,X1,LAM1,MU1,G,L1,1200);
[alpha_tikh1,ireg_corner1,kappa1] = l_curve_corner(rho1,eta1,reg_param1);
rho_corner1=rho1(ireg_corner1);
eta_corner1=eta1(ireg_corner1);

disp(['1st-order reg corner is:  ',num2str(alpha_tikh1)]);

% get the desired model
m1=m1s(:,ireg_corner1);

% plot first-order L curve
figure(3)
clf
loglog(rho1,eta1,'k-');
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Seminorm ||Lm||_2')
bookfonts

% plot the corner on the L-curve
hold on
H=loglog(rho_corner1,eta_corner1,'ko');
set(H,'markersize',16)
hold off
axis tight
print -deps2 c4flcurve1.eps

disp('Displaying the L curve for first-order regularization (fig. 3)')

% plot the first-order recovered model and the true model
figure(4)
clf
plotconst(m1*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'k--')
hold off
xlabel('Depth (m)');
ylabel('Slowness (s/km)');
bookfonts
print -deps2 c4fmtik1.eps

disp('Displaying the recovered model using first-order regularization (fig. 4)')

% apply second-order Tikhonov regularization
L2 = get_l_rough(N,2);
[U2,V2,X2,LAM2,MU2]=gsvd(G,L2);

% apply the L curve criteria to the second-order regularization problem
[rho2,eta2,reg_param2,m2s]=l_curve_tikh_gsvd(U2,t,X2,LAM2,MU2,G,L2,1200);
%second-order corner from visual inspection (l_curve_corner routine fails
%for this example)
%[alpha_tikh2,ireg_corner2,kappa2] = l_curve_corner(rho2,eta2,reg_param2);
%ireg_corner2=955;
ireg_corner2=890;
alpha_tikh2=reg_param2(ireg_corner2);
rho_corner2=rho2(ireg_corner2);
eta_corner2=eta2(ireg_corner2);

disp(['2nd-order reg corner is:  ',num2str(alpha_tikh2)]);

% get the desired model
m2=m2s(:,ireg_corner2);

%plot second-order L curve
figure(5)
clf
loglog(rho2,eta2,'k-');
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Seminorm ||Lm||_2')
bookfonts

%plot the corner on the L-curve
hold on
H=loglog(rho_corner2,eta_corner2,'ko');
set(H,'markersize',16)
hold off
axis tight
print -deps2 c4flcurve2.eps

disp('Displaying the L curve for second-order regularization (fig. 5)')

% plot the second-order recovered model and the true model
figure(6)
clf
plotconst(m2*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'k--')
hold off
xlabel('Depth (m)')
ylabel('Slowness (s/km)')
bookfonts
print -deps2 c4fmtik2.eps

disp('Displaying the recovered model using second-order regularization (fig. 6)')

% get the filter factors
gamma1=sqrt(diag(LAM1'*LAM1))./sqrt(diag(MU1'*MU1));
gamma2=sqrt(diag(LAM2'*LAM2))./sqrt(diag(MU2'*MU2));

% plot the filter factors
figure(7)
clf
semilogy(1:length(gamma1),gamma1.^2./(gamma1.^2+alpha_tikh1^2),'ko');
hold on
semilogy(1:length(gamma2),gamma2.^2./(gamma2.^2+alpha_tikh2^2),'k*');
hold off
H=text(15,0.02,'1^{st} order');
set(H,'FontSize',18);
H=text(28,0.00001,'2^{nd} order');
set(H,'FontSize',18);
xlabel('index, i');
ylabel('Filter Factors');
bookfonts
print -deps2 c4tikhfilt.eps

disp(['Displaying the filter factors for first and second-order'...
    ' regularization (fig. 7)'])

% display the recovered model misfits
dinterp=dd/2:dd:1000-dd/2;
strueresamp=slowfunctionsmooth(dinterp)';

disp(['1st-order model 2-norm misfit is ',...
    num2str(norm(m1-strueresamp)),' s/km']);
disp(['2nd-order model 2-norm misfit is ',...
    num2str(norm(m2-strueresamp)),' s/km']);

% Example 4.5
% examine resolution
%
% compute the resolution matrix for first-order
f1=gamma1.^2./(gamma1.^2+alpha_tikh1^2);
for i=1:length(f1);
    if isnan(f1(i))
        f1(i)=1;
    end
end
        
F1=diag([f1;ones(N-length(f1),1)]);
R1=inv(X1)'*F1*X1';

% compute the resolution matrix for second-order
f2=gamma2.^2./(gamma2.^2+alpha_tikh2^2);
for i=1:length(f2);
    if isnan(f2(i))
        f2(i)=1;
    end
end
F2=diag([f2;ones(N-length(f2),1)]);
R2=inv(X2)'*F2*X2';

% spike resolution for example figure
spike=zeros(N,1);
spike(N/2)=1;

% the models that would be recovered
r1=R1*spike;
r2=R2*spike;

% plot the recovered spike models
figure(8)
clf
subplot(2,1,1)
plotconst(r1,0,maxdepth);
ylim([0 0.72])
xlim([0 1000])
bookfonts
H=text(700,.15,'1^{st} Order');
set(H,'FontSize',18);
ylim([0 .2])

subplot(2,1,2)
plotconst(r2,0,maxdepth);
ylim([0 0.72])
xlim([0 1000])
xlabel('Depth (m)')
bookfonts
H=text(700,.15,'2^{nd} Order');
set(H,'FontSize',18);
ylim([0 .2])
print -deps2 c4fspikeres.eps

disp('Displaying recovered spike models (fig. 8)')
