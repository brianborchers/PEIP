% Example 7.4
% from Parameter Estimation and Inverse Problems, 3nd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% This script generates the G matrix, mtrue, dtrue, and noisy d 
% then solves the inverse problem for a smooth true model and for a second
% true model that has rapid changes

% make sure we have a clean environment
clear 
rand('state',0);
randn('state',0);

% set up input for getrep.m 
xsample=[0.01,25.05,50:10:250,275,300];
tsample=300;

% params is fed to kernel.m with specific x and t parameters
% params = [velocity, diffusion]
params=[1.0, 1.0];

% the time the model will occupy
nt=100;
tmax=250.;
tmin=0.01;
% the discrete time points
t=linspace(tmin,tmax,nt);

% get true source concentration as a function of t (smooth)
cin=cintrue(t);

% Make the G matrix.
G=getrep('kernel',t,xsample,tsample,params);

% generate sampled data
cexact=G*cin';
csample=cexact+0.001*randn(size(cexact));

% rename some things so the rest of this is in a form we're used to
mtrue=cin';
d=csample;
dtrue=cexact;

% Plot out the true model
figure(1)
clf
plot(t,cin,'k-');
xlabel('Time (days)');
ylabel('Concentration');
bookfonts
disp('Displaying the true model (fig. 1)')
% print -deps2 c7fcintrue.eps

% Plot the downstream samples
figure(2)
clf
plot(xsample,d,'k-');
xlabel('Distance');
ylabel('Concentration');
bookfonts
disp('Displaying the simulated data (fig. 2)')
print -deps2 c7fcout.eps

% Get the problem size.
[m,n]=size(G);

% Use TV regularization to solve the problem (smooth true model)
% Setup the L matrix.
L=get_l_rough(n,1);

% The alphas to compute the curve at
alphas=10.^(-8:.25:-2)';

% Compute the tradeoff curve.
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
for i=1:length(alphas)
  mreg(:,i)=irlsl1reg(G,d,L,alphas(i));
  rho(i)=norm(G*mreg(:,i)-d);
  eta(i)=norm(L*mreg(:,i),1);
end

% get the curve corner
[reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,alphas);

% Plot the tradeoff curve.
figure(3)
clf
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Seminorm ||Lm||_{1}');
% mark and label the corner
hold on
loglog(rho(ireg_corner),eta(ireg_corner),'ko','MarkerSize',18)
H=text(rho(ireg_corner),eta(ireg_corner),...
    ['   ',num2str(alphas(ireg_corner),'%5.1e')]);
set(H,'FontSize',18);
% label to each side of the corner
H=text(rho(20),eta(20),['   ',num2str(alphas(20),'%5.1e')]);
set(H,'Fontsize',18);
H=text(rho(5),eta(5),['      ',num2str(alphas(5),'%5.1e')]);
set(H,'Fontsize',18);
hold off
bookfonts

% axis([20 200 1 1e4])

disp('Displaying the total variation L-curve (fig. 3)')
print -depsc2 c7contamtvlcurve.eps

% Plot the suite of solutions
figure(4)
clf
hold on
for i=2:2:length(alphas)
  plot(t,mreg(:,i)*.6+log10(alphas(i)),'k');
end
plot(t,mreg(:,ireg_corner)*.6+log10(alphas(ireg_corner)),'k','LineWidth',2);
xlabel('Time (days)');
ylabel('log_{10}(\alpha)');
axis tight
bookfonts

disp('Displaying a suite of total variation recovered models (fig. 4)')
print -depsc2 c7contamtvsols.eps

% Plot the TV solution
figure(5)
clf
plot(t,mreg(:,ireg_corner),'k',t,mtrue,'k--');
xlabel('Time (days)');
ylabel('m(t)');
bookfonts

disp('Displaying the L-curve corner total variation recovered model (fig. 5)')
print -depsc2 c7contamtvsol.eps

% get a second true source concentration as a function of t (one that
% actually has rapid changes
cin=cintrue_TV(t);

% Generate sampled data
cexact=G*cin';
csample=cexact+0.001*randn(size(cexact));

% rename some things so the rest of this is in a form we're used to
mtrue=cin';
d=csample;
dtrue=cexact;

% Plot out the true model
figure(6)
clf
plot(t,cin,'k-');
xlabel('Time (days)');
ylabel('Concentration');
axis([0 250 0 0.8]);
bookfonts

disp('Displaying the true model (fig. 6)')
print -deps2 c7fcintrue_TV.eps

% Plot the downstream samples
figure(7)
clf
plot(xsample,d,'k-');
xlabel('Distance');
ylabel('Concentration');
bookfonts

disp('Displaying the simulated data (fig. 7)')
print -deps2 c7fcout2.eps

% 
%Solve the problem again, using TV regularization (true model with rapid
%changes)
% Setup the L matrix.
L=get_l_rough(n,1);

% The alphas to compute the curve at
alphas=10.^(-8:.25:-2)';

% Compute the tradeoff curve.
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
for i=1:length(alphas)
  mreg(:,i)=irlsl1reg(G,d,L,alphas(i));
  rho(i)=norm(G*mreg(:,i)-d);
  eta(i)=norm(L*mreg(:,i),1);
end

% get the curve corner
[reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,alphas);

% Plot the tradeoff curve.
figure(8)
clf
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Seminorm ||Lm||_{1}');
% mark and label the corner
hold on
loglog(rho(ireg_corner),eta(ireg_corner),'ko','MarkerSize',18)
H=text(rho(ireg_corner),eta(ireg_corner),...
    ['   ',num2str(alphas(ireg_corner),'%5.1e')]);
set(H,'FontSize',18);
% label to each side of the corner
H=text(rho(20),eta(20),['   ',num2str(alphas(20),'%5.1e')]);
set(H,'Fontsize',18);
H=text(rho(5),eta(5),['      ',num2str(alphas(5),'%5.1e')]);
set(H,'Fontsize',18);
hold off
bookfonts

disp('Displaying the total variation L-curve (fig.8)')
print -depsc2 c7contamtvlcurve2.eps

% Plot the suite of solutions
figure(9)
clf
hold on
for i=2:2:length(alphas)
  plot(t,mreg(:,i)*.6+log10(alphas(i)),'k');
end
plot(t,mreg(:,ireg_corner)*.6+log10(alphas(ireg_corner)),'k','LineWidth',2);
xlabel('Time (days)');
ylabel('log_{10}(\alpha)');
axis tight
bookfonts

disp('Displaying a suite of total variation recovered models (fig. 9)')
print -depsc2 c7contamtvsols2.eps

% Plot the denoised solution
figure(10)
clf
plot(t,mreg(:,ireg_corner),'k',t,mtrue,'k--');
xlabel('Time (days)');
ylabel('m(t)');
axis([0 250 0 0.8]);
bookfonts

disp('Displaying the L-curve corner total variation recovered model (fig. 10)')
print -depsc2 c7contamtvsol2.eps
