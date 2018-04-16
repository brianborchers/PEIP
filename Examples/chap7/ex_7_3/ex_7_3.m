% Example 7.3
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% This script demonstrates compressive sensing by denoising a
% 1-d signal that happens to be sparse with respect to the discrete cosine transform basis.

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% We'll use 100 Hz sampling, and a ten second time period to get 1001
% signal points in the series.
r=100;
deltat=1/r;
t=0.0:deltat:10;
n=length(t);

% For our input signal, we'll use a mixture of 2 sine waves.  To
% avoid spectral leakage, we'll multiply this by a Hanning window.
f1=25;               % first sine wave at 25 Hz.
f2=35;               % second sine wave at 35 Hz.
%construct the signal
mtrue=(5*cos(2*pi*f1*t)+2*cos(2*pi*f2*t))';
%taper the signal
mtrue=0.5*(1 - cos(2*pi*(1:n)'/(n+1))).*mtrue;
xtrue=dcost(mtrue);

% Plot the input signal.
figure(1)
clf
plot(t,mtrue,'k');
xlabel('Time (s)');
ylabel('m(t)');
bookfonts

disp('Displaying the true input signal (fig. 1)')
print -depsc2 c7fcssignal.eps

% We'll take m=100 measurements.
m=100;

% Take m observations from the n data points in the signal by
% using a G matrix that contains random aggregations of the data points
G=randn(m,n);

% The data consists of only those m samples from the original signal.
dtrue=G*mtrue;

% Add noise to the data.
noise=5;
d=dtrue+noise*randn(size(dtrue));

% First, use Tikhonov regularization to denoise the signal directly.
L=get_l_rough(n,2);
alphas=10.^(-12:5)';

% Compute and plot the l-curve.
rhoT=zeros(length(alphas),1);
etaT=zeros(length(alphas),1);
mreg=zeros(n,length(alphas));

% disable printing of expected warnings
warning('off', 'MATLAB:rankDeficientMatrix');

% Loop through the alphas finding regularized solutions.
for i=1:length(alphas)
  mreg(:,i)=[G; alphas(i)*L]\[d; zeros(n-2,1)];
  rhoT(i)=norm(G*mreg(:,i)-d);
  etaT(i)=norm(L*mreg(:,i));
end

% reenable printing of expected warnings
warning('off', 'MATLAB:rankDeficientMatrix');

% Based on magnitude of the noise, which should be sigma*sqrt(n)~150, we pick: 
ireg_cornerT=14;

figure(2)
clf
semilogx(rhoT,etaT,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Seminorm ||Lm||_{2}');
bookfonts

% Show the point for the discrepancy principle-determined value of alpha.
hold on
H=plot(rhoT(ireg_cornerT),etaT(ireg_cornerT),'ko');
set(H,'markersize',16)

disp('Displaying the L-curve for second order regularization (fig. 2) and its discrepancy principal value')

% Plot the Tikhonov regularized solution selected from the discrepancy principle
figure(3)
clf
plot(t,mreg(:,ireg_cornerT),'k');
xlabel('Time (s)');
ylabel('m(t)')
ylim([-10 10]);
bookfonts

disp(['Displaying the discrepancy principle determined regularized'...
    ' model (fig. 3)'])
print -depsc2 c7fcstiksol.eps

% Plot a suite of the Tikhonov regularized models
figure(4)
clf
hold on
for i=1:length(alphas)
  plot(t,0.1*mreg(:,i)+log10(alphas(i)),'k');
end
plot(t,0.1*mreg(:,ireg_cornerT)+log10(alphas(ireg_cornerT)),'k','LineWidth',2);
xlabel('Time (s)');
ylabel('log_{10}(\alpha)');
bookfonts
axis tight

disp('Displaying a suite of regularized models (fig. 4)')

% Now, use the sparsity regularization to get the signal back.  

% We'll need an n by n identity matrix.
I=eye(n);

% First construct an explicit basis for the discrete cosine transform.
% where the columns of W are the discrete cosine transforms of the standard basis
W=dcost(I);

% The A matrix is given by A=G*W
A=G*W;

alphas=10.^(-12:5)';

% disable printing of expected warnings
warning('off', 'MATLAB:nearlySingularMatrix');

% Compute the tradeoff curve.
xcs=zeros(n,length(alphas));
rhocs=zeros(length(alphas),1);
etacs=zeros(length(alphas),1);
for i=1:length(alphas)
  xcs(:,i)=irlsl1reg(A,d,I,alphas(i));
  rhocs(i)=norm(A*xcs(:,i)-d);
  etacs(i)=norm(xcs(:,i),1);
end

% reenable printing of expected warnings
warning('on', 'MATLAB:nearlySingularMatrix');

% Picking the corner by hand from the unusual L curve
ireg_cornercs=2;

% Plot the sparse tradeoff curve.
figure(5)
clf
loglog(rhocs,etacs,'k-');
xlabel('Residual Norm ||GWx - d||_{2}');
ylabel('Solution Norm ||x||_{1}');
bookfonts

disp('Displaying the sparse tradeoff curve (fig. 5)')
% Show the point for the hand picked alpha value
hold on
H=plot(rhocs(ireg_cornercs),etacs(ireg_cornercs),'ko');
set(H,'markersize',16)

disp('Displaying the trade off curve for the sparse solution (fig. 5)')

% reconstruct and plot the regularized solutions.
for i=1:length(alphas)
  mcs(:,i)=W*xcs(:,i);
end

% plot the selected regularized solution
figure(6)
clf
plot(t,mcs(:,ireg_cornercs),'k');
xlabel('Time (s)');
ylabel('m(t)');
bookfonts
disp('Displaying the selected sparse solution (fig. 6)')
print -depsc2 c7fcssol.eps

% plot the suite of regularized solutions
figure(7)
clf
hold on
for i=1:2:length(alphas)
  plot(t,mcs(:,i)/max(mcs(:,i))+log10(alphas(i)),'k');
end
plot(t,...
    mcs(:,ireg_cornercs)/max(mcs(:,ireg_cornercs))+...
    log10(alphas(ireg_cornercs)),'k','LineWidth',2);
axis tight
xlabel('Time (s)');
ylabel('log_{10}(\alpha)');
bookfonts

disp('Displaying a suite of regularized sparse solutions (fig. 7)')
