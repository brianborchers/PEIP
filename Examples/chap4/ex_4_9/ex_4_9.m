% Example 4.9
% from Parameter Estimation and Inverse Problems, 3nd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% This script generates the G matrix, mtrue, dtrue, and noisy d 
% then solves the inverse problem in various ways

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

% get true source concentration as a function of t
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
xlabel('Time');
ylabel('Concentration');
bookfonts

disp('Displaying the true model (fig. 1)')
print -deps2 c4fcintrue.eps

% Plot the downstream samples
figure(2)
clf
plot(xsample,d,'k-');
xlabel('Distance');
ylabel('Concentration');
bookfonts

disp('Displaying the simulated data (fig. 2)')
print -deps2 c4fcout.eps

% Get the problem size.
[m,n]=size(G);

% First, get the least squares solution.
mls=pinv(G)*d;

% plot the least squares solution
figure(3)
clf
plot(t,mls,'k-',t,cin,'k--');
ylabel('Concentration');
xlabel('Time');
bookfonts

disp('Displaying least squares recovered model (fig. 3)')
print -deps2 c4fmls.eps

% Next, solve the problem using non-negative least squares
mlsn=lsqnonneg(G,d);


% plot the non-negative least squares solution
figure(4)
clf
plot(t,mlsn,'k-',t,cin,'k--');
xlabel('Time');
ylabel('Concentration');
bookfonts

disp('Displaying non-negative least squares recovered model (fig. 4)')
print -deps2 c4fmlsn.eps

% Next, solve the problem using BVLS
% set up the bounds
l=zeros(n,1);
u=1.1*ones(n,1);
mbvls=bvls(G,d,l,u);

% plot the bounded least squares solution
figure(5)
clf
plot(t,mbvls,'k-',t,cin,'k--');
xlabel('Time');
bookfonts
ylabel('Concentration');
disp('Displaying bounded least squares recovered model (fig. 5)')
print -deps2 c4fmbvls.eps


% Next, solve the problem using 2nd order Tikhonov regularization with BVLS
alphas=10.^(-5:0.25:1)';
rhos=zeros(size(alphas))';
eta=zeros(size(alphas))';

% Setup the L matrix.
L=zeros(n-2,n);
for i=1:n-2
  L(i,i)=1;
  L(i,i+1)=-2;
  L(i,i+2)=1;
end

% Add 0's to the end of d to get the right hand side of the regularized 
% least squares problems.
dpad=[d; zeros(n-2,1)];

% Now, compute the regularized solutions.
for i=1:length(alphas)
  A=[G; alphas(i)*L];
  m=bvls(A,dpad,l,u);
  rho(i)=norm(G*m-d);
  eta(i)=norm(L*m);
end  

% Find the corner of the L-curve.
[reg_corner,ireg_corner]=l_curve_corner(rho',eta',alphas);
disp(['The corner of the L-curve occurs at alpha=',...
        num2str(alphas(ireg_corner))]);


% Plot the L-curve.
figure(6)
clf
% plot the curve
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_2');
ylabel('Solution Seminorm ||Lm||_2');
hold on
% mark and label the corner
plot(rho(ireg_corner),eta(ireg_corner),'ko','MarkerSize',18);
H=text(rho(ireg_corner),eta(ireg_corner),...
    ['  ',num2str(alphas(ireg_corner),'%5.1e')]);
set(H,'FontSize',18);
% label to each side of the corner
H=text(rho(ireg_corner-5),eta(ireg_corner-5),...
    ['  ',num2str(alphas(ireg_corner-5),'%5.1e')]);
set(H,'FontSize',18);
H=text(rho(ireg_corner+5),1.1*eta(ireg_corner+5),...
    ['    ',num2str(alphas(ireg_corner+5),'%5.1e')]);
set(H,'FontSize',18);
bookfonts

disp('Displaying the second order L-curve (fig. 6)')
print -deps2 c4flcurve_shist.eps


% Get the corner model
alpha=alphas(ireg_corner);
A=[G; alpha*L];
mtik=bvls(A,dpad,l,u);

% plot the Tikhonov model
figure(7)
clf
plot(t,mtik,'k-',t,cin,'k--');
xlabel('Time');
ylabel('Concentration');
bookfonts

disp('Displaying the second order Tikhonov recovered model (fig. 7)')
print -deps2 c4fmtik_shist.eps


% Now, construct bounds on the minimum/maximum average value in the
% interval from t=125 to t=150.
c=zeros(n,1);
c(51:60)=0.1*ones(10,1);
[xmin,xmax]=blf2(G,d,c,0.01,l,u);

disp(['Average concenteration between ', num2str(c'*xmin), ' and ',...
    num2str(c'*xmax),'.'])

