% Example 7.1
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% This script demonstrates deconvolution of a spike train by
% minimization of the 1-norm of the model subject to a quadratic
% constraint on data fit.  

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% We'll use 100 Hz sampling, and a ten second time period.
r=100;
deltat=1/r;
t00=0.5;
t=(0:999)*deltat-0.5;
noise=0.002;
m=zeros(length(t),1);

%
% Test model -- response of a single reflector with multiples every 1.5 s
% that flip in sign.

%pre-source time
T0=t00*r;
%two-way travel time in samples
tau=150;
%number of multiples
nmult=5;
%initial amplitude
A=10;
%reflection coefficient
R=0.4;

%source spike
m(T0)=A;
for i=1:nmult
    m(i*tau+T0)=A*R^i*(-1)^i;
end

% Plot the true model
figure(1)
clf
plot(t,m,'k');
xlabel('Time (s)');
ylabel('Impulse Response Amplitude');
bookfonts
xlim([-0.5 10])

disp('Displaying the true model (fig. 1)')
print -deps2 c7fpulsetarget.eps 

% For our impulse response we'll use g(t)=exp(-5t)*sin(10t)
tshort=0.0:deltat:2.0;
g=exp(-5*tshort).*sin(10.*tshort);

% Plot the impulse response.
figure(2)
clf
plot(tshort,g,'k');
axis tight
xlabel('Time (s)');
ylabel('Source Amplitude');
bookfonts

disp('Displaying the instrument impulse response (fig. 2)')

% Perform the convolution to calculate the noise free data
d=conv(g,m)*deltat;

% Cut the data to 10 seconds.
N=1000;
d=d(1:N);
% Also get the noisy data vector
dn=d+noise*randn(size(d));


% Plot the clean data.
figure(3)
clf
plot(t,d,'k');
xlim([-0.5 10])
axis tight
xlabel('Time (s)');
ylabel('Seismic Amplitude');
bookfonts

disp('Displaying the noise free data (fig. 3)')

% Plot the noisy data.
figure(4)
clf
plot(t,dn,'k');
xlim([-0.5 10])
axis tight
xlabel('Time (s)');
ylabel('Seismic Amplitude');
bookfonts

disp('Displaying the noisy data (fig. 4)')
print -deps2 c7fnoisypulses.eps


% Try standard zeroth-order Tikhonov regularization.
% Create the model
G=zeros(length(dn)+length(g),length(t));
for i=1:length(t)
  G(i:i+length(g)-1,i)=g'*deltat;
end
G=G(1:N,:);
GGT=G'*G;

L=speye(N);

alphas=10.^(-3.25:.25:5);
rhoT=zeros(length(alphas),1);
etaT=zeros(length(alphas),1);
mtik=zeros(N,length(alphas));

for i=1:length(alphas)
    mtik(:,i)=(GGT+alphas(i)^2*L)\(G'*dn);
    rhoT(i)=norm(G*mtik(:,i)-dn);
    etaT(i)=norm(mtik(:,i));
end

%get the corener of the L-curve
[~,ireg_cornerT,~]=l_curve_corner(rhoT,etaT,alphas);
rho_cornerT=rhoT(ireg_cornerT);
eta_cornerT=etaT(ireg_cornerT);


% plot the zeroth order L-curve
figure(5)
clf
plot(rhoT,etaT,'k.-');
xlabel('Residual Norm ||Gm - d||_{2}');
ylabel('Solution Norm ||m||_{2}');
bookfonts
hold on
% mark and label the corner
H=plot(rho_cornerT,eta_cornerT,'ko');
set(H,'markersize',16)
H=text(rho_cornerT+0.01,eta_cornerT+1,num2str(alphas(ireg_cornerT)));
set(H,'Fontsize',18);
% label to both sides on the L-curve
H=text(rhoT(2)+0.005,etaT(2)+1,num2str(alphas(2),3));
set(H,'Fontsize',18);
H=text(rhoT(10)+0.01,etaT(10)+2,num2str(alphas(10),3));
set(H,'Fontsize',18);
axis tight

disp('Displaying zeroth-order Tikhonov L-curve (fig. 5)')
print -deps2 c7ftikLcurve.eps

% Plot the Tikhonov zeroth-order inverse solution.
figure(6)
clf
plot(t,mtik(:,ireg_cornerT),'k');
xlim([-0.5 10])
axis tight
xlabel('Time (s)');
ylabel('Impulse Response Amplitude');
bookfonts

disp('Displaying the zeroth order Tikhonov solution (fig. 6)')
print -deps2 c7ftikpulse.eps

% Now, perform the sparse deconvolution.
disp('Note: Performing sparse deconvolution; this computation may take several minutes on most computers.')

%
% Sparsify G.
%
G=sparse(G);

%examine a logarithmically spaced range of regularization parameters
alphas=10.^(-6:.125:-2);
%
% Get the Lipschitz constant.
%
Lip=2.05*normest(G)*normest(G');
%
% Allocate space for the results.
%
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
msparse=zeros(length(t),length(alphas));
mreg=zeros(length(t),1);
%
% Loop over the values of alpha.
%
for i=1:length(alphas);
  fprintf('Starting with alpha=%e\n',alphas(i));
  mreg=fista(G,dn,alphas(i),Lip,[],1.0e-5,mreg);
  msparse(:,i)=mreg;
  rho(i)=norm(G*msparse(:,i)-dn);
  eta(i)=norm(msparse(:,i),1);
end

% Find the corner of the Tikhonov L-curve
[reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,alphas);

% Plot the sparse deconvolution L-curve.
figure(7)
clf
loglog(rho,eta,'k.-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||m||_{1}');
bookfonts
hold on
% mark and label the corner
H=loglog(rho(ireg_corner),eta(ireg_corner),'ko');
set(H,'markersize',8)
H=text(rho(ireg_corner),1.1*eta(ireg_corner),...
    ['    ',num2str(alphas(ireg_corner),'%5.1e')]);
set(H,'Fontsize',18);
% label to both sides on the L-curve
H=text(rho(6),eta(6),['  ',num2str(alphas(6),'%5.1e')]);
set(H,'Fontsize',18);
H=text(rho(28),eta(28),['    ',num2str(alphas(28),'%5.1e')]);
set(H,'Fontsize',18);
axis([3e-2 .3 1 1e3])

disp('Displaying the 1-norm zeroth-order L-curve (fig. 7)')
print -deps2 c7fsparsepulseL.eps

% plot the solution for the corner
figure(8)
clf
plot(t,msparse(:,ireg_corner),'k');
xlabel('Time (s)');
ylabel('Impulse Response Amplitude');
bookfonts
axis tight

disp('Displaying 1-norm regularized solution (fig. 8)')
print -deps2 c7fsparsepulse.eps
%
% Output some information about the solution.
%
fprintf('Corner alpha=%e\n',alphas(ireg_corner));
fprintf('one-norm of m is %e\n',norm(msparse(:,ireg_corner),1));
fprintf('|| Gm-d || is %e\n',norm(G*msparse(:,ireg_corner),2));
fprintf('Weighted objective is %e\n',norm(G*msparse(:,ireg_corner),2)+...
	alphas(ireg_corner)*norm(msparse(:,ireg_corner),1));

% Plot a selection of the sparse solutions sorted by regularization
% parameter
figure(9)
clf
hold on
for i=3:2:length(alphas)
  plot(t,msparse(:,i)*.06+log10(alphas(i)),'k');
  xlabel('Time (s)');
  ylabel('log_{10} \alpha');
  bookfonts
end
plot(t,msparse(:,ireg_corner)*.06+log10(alphas(ireg_corner)),'k','LineWidth',2);
hold off
axis tight

disp('Displaying a range of 1-norm regularized solutions (fig. 9)')
print -deps2 c7fsparsepulsesols.eps

%examine solutions in an animation
disp(['Animating the 1-norm regularized models recovered'...
' with increasing alpha (fig. 10)'])
figure(10)
clf
for i=1:length(alphas)
  plot(t,msparse(:,i),'k');
  axis tight
  xlabel('Time (s)');
  ylabel('Impulse Response Amplitude');
  title(['log_{10}(\alpha): ',num2str(log10(alphas(i)))]);
  bookfonts
  pause(0.1)
end
% 
% Output information about the corner solution.
%
alpha=alphas(ireg_corner);
mreg=msparse(:,ireg_corner);
obj=norm(G*mreg-dn,2)^2+alpha*norm(mreg,1);
fprintf('For corner solution:\n');
fprintf('alpha=%e\n',alpha);
fprintf('Weighted objective value is %e\n',obj);
fprintf('rho=%e, eta=%e\n',[rho(ireg_corner); eta(ireg_corner)]);
