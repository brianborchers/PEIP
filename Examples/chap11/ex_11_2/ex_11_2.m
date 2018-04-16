% Example 11.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% Reinitialize everything.
%
clear
rand('seed',0);
randn('seed',0);
%
% Get the problem.
%
load shawexamp.mat
%
% Set the data covariance matrix.
%
CD=1.0e-12*eye(20);
%
% Setup the prior.
%
mprior=0.5*ones(20,1);
CM=0.25*eye(20); 
%
% Compute the posterior distribution.
%
[covmp,mmap]=bayes(G,mprior,CM,dn,CD);
%
% Compute a second solution with a more diffuse prior.
%
mprior2=mprior;
CM2=100*CM;
[covmp2,mmap2]=bayes(G,mprior2,CM2,dn,CD);
%
% Now, produce some plots.
%
figure(1)
clf
plotconst(spike,-pi/2,pi/2);
axis([-2 2 -0.5 1.5]);
hold on
plotconstc(mmap,-pi/2,pi/2,'k--');
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('Target Model','MAP Solution');
print -deps2 c11fmmap.eps
disp('Displaying the true model and a MAP solution based on an uninformative prior (fig. 1)');
%
% Now, with error bars.
%
figure(2)
clf
plotconst(mmap,-pi/2,pi/2);
axis([-2 2 -1.0 1.5]);
hold on
plotconstc(mmap+1.96*sqrt(diag(covmp)),-pi/2,pi/2,'k--');
plotconstc(mmap-1.96*sqrt(diag(covmp)),-pi/2,pi/2,'k--');
hold off
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('MAP Solution','95% Prob Interval');
print -deps2 c11fmmapb.eps
disp('Displaying the MAP solution with 95% confidence interval (fig. 2)');
%
% Generate a random solution.
%
figure(3)
clf
%
% Use the library function simmvn to generate a random solution.
%
mmapsims=simmvn(mmap,covmp);
plotconst(mmapsims,-pi/2,pi/2);
xlabel('\theta');
ylabel('Intensity');
bookfonts
print -deps2 c11fmmapsims.eps
disp('Displaying randomly selected realization (fig. 3)');
%
% Now, with error bars for the more diffuse prior.
%
figure(4)
clf
plotconst(mmap2,-pi/2,pi/2);
hold on
plotconstc(mmap2+1.96*sqrt(diag(covmp2)),-pi/2,pi/2,'k--');
plotconstc(mmap2-1.96*sqrt(diag(covmp2)),-pi/2,pi/2,'k--');
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('MAP Solution','95% Prob Interval');
print -deps2 c11fmmap2.eps
disp('Displaying the true model and a MAP solution based on an even less uninformative prior (fig. 4)');

%
% Setup a more restrictive prior.
%
mprior=zeros(20,1);
mprior(6:14)=0.5*(1 - cos(2*pi*(1:9)'/(10)));
 
CM=50*diag(0.00125*ones(20,1).*(0.5*(1 - cos(2*pi*(1:20)/(21))))'.^2); 

%
% Compute the posterior distribution.
%
[covmp,mmap]=bayes(G,mprior,CM,dn,CD);
%
figure(5)
clf
plotconst(mprior,-pi/2,pi/2);
hold on
plotconstc(mprior+1.96*sqrt(diag(CM)),-pi/2,pi/2,'k--');
plotconstc(mprior-1.96*sqrt(diag(CM)),-pi/2,pi/2,'k--');
hold off
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('Prior Distribution');
print -deps2 c11fpriorr.eps
disp('Displaying the prior distribution and a 95% confidence interval around it (fig. 5)');

figure(6)
clf
plotconst(spike,-pi/2,pi/2);
axis([-2 2 -0.5 1.5]);
hold on
plotconstc(mmap,-pi/2,pi/2,'k--');
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('Target Model','MAP Solution');
print -deps2 c11fmmapr.eps
disp('Displaying the true model and a MAP solution based on prior that expects a central spike (fig. 6)');
%
% Now, with error bars.
%
figure(7)
clf
plotconst(mmap,-pi/2,pi/2);
axis([-2 2 -1.0 1.5]);
hold on
plotconstc(mmap+1.96*sqrt(diag(covmp)),-pi/2,pi/2,'k--');
plotconstc(mmap-1.96*sqrt(diag(covmp)),-pi/2,pi/2,'k--');
hold off
xlabel('\theta');
ylabel('Intensity');
bookfonts
legend('MAP Solution','95% Prob Interval');
print -deps2 c11fmmapbr.eps
disp('Displaying the MAP solution with 95% confidence interval (fig. 7)');
%
% Generate a random solution.
%
figure(8)
mmapsims=simmvn(mmap,covmp);
plotconst(mmapsims,-pi/2,pi/2);
bookfonts
xlabel('\theta');
ylabel('Intensity');
bookfonts
disp('Displaying a random realization of the restricted model (fig. 8)');
