% Examples 4.1, 4.2, 4.3, and 4.8
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('seed',0);
randn('seed',0);

% Load in the Shaw problem for m=n=20;
load shawexamp.mat

% Example 4.1
% Calculate the svd
[U,S,V]=svd(G);

% First, calculate and plot the L-curve, and find its corner.
s=diag(S);

[rho,eta,reg_param]=l_curve_tikh_svd(U,s,dspiken,1000);

% estimate the L-curve corner in log-log space using a spline fit to find where
% maximum negative curvature occurs(requires Spline Toolbox)
[alpha_tikh,ireg_corner] = l_curve_corner(rho,eta,reg_param);
rho_corner=rho(ireg_corner);
eta_corner=eta(ireg_corner);

disp(['alpha from L curve corner is ', num2str(alpha_tikh)]);

% get the spike solution corresponding to the L-curve corner
m_tikh=(G'*G+alpha_tikh^2*eye(20))\G'*dspiken;

% get and display the residual using the L-curve solution
r_spike=norm(G*m_tikh-dspiken);

disp(['Residual norm for L-curve solution using Tikhonov regularization: ',num2str(r_spike)]);

% plot L curve and add the corner marker
figure(1)
clf
loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Norm ||m||_2')
bookfonts
axis tight
hold on
H=loglog(rho_corner,eta_corner,'ko');
set(H,'markersize',16)
hold off
print -deps2 c4flcurve0

disp('Displaying the L-curve (fig. 1)')

% plot the L curve predicted model
figure(2)
clf
plotconst(m_tikh,-pi/2,pi/2);
xlabel('\theta');
ylabel('Intensity');
bookfonts
ylim([-.2 0.5])
print -deps2 c4fmtik.eps

disp('Displaying the L-curve predicted model (fig. 2)')


% Use the discrepancy principle to get a second solution.
% find the regularization value, alpha_disc, for rho=discrep by interpolation 
% of the L-curve
discrep=(1e-6)*sqrt(20);
alpha_disc=interp1(rho,reg_param,discrep);

disp(['alpha from the discrepancy principle is ', num2str(alpha_disc)])

% get the model and residual
m_disc=(G'*G+alpha_disc^2*eye(20))\G'*dspiken;
r_spike_disc=norm(G*m_disc-dspiken);

disp(['Residual norm for discrepancy principle solution using Tikhonov',...
    ' regularization: ',num2str(r_spike_disc)]);

% plot the discrepancy principle predicted model
figure(3)
clf
plotconst(m_disc,-pi/2,pi/2);
xlabel('\theta');
ylabel('Intensity');
bookfonts
ylim([-.2 0.5])
print -deps2 c4fmdisc.eps

disp('Displaying the discrepancy principle predicted model (fig. 3)')


% get the values for the Picard plot
[utd,utd_norm]=picard_vals(U,s,dspiken);

% Produce the Picard plot.
figure(4)
clf
x_ind=1:length(s);
semilogy(x_ind,s,'k-',x_ind,abs(utd),'k.',x_ind,abs(utd_norm),'ko')
legend('s_i','|u_i^Td|','|u_i^Td|/s_i','Location','SouthWest');
xlabel('i')
bookfonts
axis tight
print -deps2 c4fpicard.eps

disp('Displaying the Picard plot (fig. 4)')

% Example 4.2
% Now, examine the resolution using a noise-free spike test for alpha_tikh
rdspike=(G'*G+alpha_tikh^2*eye(20))\G'*dspike;

% plot the noise-free spike model for the L curve
figure(5)
clf
plotconst(rdspike,-pi/2,pi/2);
ylim([-.2 0.5])
xlabel('\theta');
ylabel('Intensity');
bookfonts

disp(['Displaying the predicted model from the noise free data using the'...
    ' L-curve criteria (fig. 5)'])

% Now, examine the resolution using a noise-free spike test for alpha_disc
rdspike=(G'*G+alpha_disc^2*eye(20))\G'*dspike;

% plot the noise-free spike model for the discrepancy principle
figure(6)
clf
plotconst(rdspike,-pi/2,pi/2);
ylim([-.2 0.5])
xlabel('\theta');
ylabel('Intensity');
bookfonts
print -deps2 c4fmdisc_noise_free.eps

disp(['Displaying the predicted model from the noise free data using the'...
    ' discrepancy principle (fig. 6)'])

% Compute the resolution matrix for alpha_disc.
Ginv_disc=(G'*G+alpha_disc^2*eye(20))\G';
R_disc=Ginv_disc*G;

disp('diagonal resolution elements:')
diag(R_disc)

% Plot the resolution matrix
figure(7)
clf
imagesc(R_disc)
axis square
colormap('gray')
colorbar
xlabel('j');
ylabel('i')
bookfonts
print -deps2 c4shaw_res.eps

disp(['Displaying the resolution matrix for the discrepancy principle'...
    ' (fig. 7)'])

% Example 4.3
% get the covariance of the discrepancy principle solution
covm_disc=Ginv_disc*1.0e-6^2*eye(20)*Ginv_disc';

% Now, produce a plot of the discrepancy principle solution with "error
% bars" versus reality.  
figure(8)
clf
H1=plotconstc(spike,-pi/2,pi/2,'k-');
hold on
H2=plotconstc(m_disc,-pi/2,pi/2,'k--');
H3=plotconstc(m_disc+1.96*sqrt(diag(covm_disc)),-pi/2,pi/2,'k:');
legend(H1,'m_{true}','m_{disc}','95%');
H4=plotconstc(m_disc-1.96*sqrt(diag(covm_disc)),-pi/2,pi/2,'k:');
xlabel('\theta');
ylabel('Intensity');
bookfonts
print -deps2 c4fbias.eps

disp(['Displaying the discrepancy principle predicted model with error'...
    ' bars (fig. 8)'])

% Example 4.8
% Tikhonov's theorem calculations
w=(G')\spike;
disp(['norm(G^T\d_spike) is ', num2str(norm(w))]);


% This is hopelessly large, so we'll give up and try a smoother model.
% generate the true model, noise free and noisy data
w=[ones(3,1); zeros(17,1)];
smoothmod=G'*w;
smoothmodd=G*smoothmod;
smoothmoddn=smoothmodd+1.0e-6*randn(20,1);

% Choose alpha (alphahat) based on theorem 4.2
p=1;
gamma=1/2;
delta=sqrt(20)*1.0e-6;
Delta=delta/norm(w);
alphahat=(Delta/2*gamma*1)^(1/(1+p));

% the predicted smooth model with the above alpha
malphahat=(G'*G+alphahat^2*eye(20))\G'*smoothmoddn;

disp(['norm(malphahat-smoothmod) is ', num2str(norm(malphahat-smoothmod))]);

% plot the true smooth model
figure(9)
clf
plotconst(smoothmod,-pi/2,pi/2);
xlabel('\theta');
ylabel('Intensity');
bookfonts
print -deps c4fsmoothmod.eps

disp(['Displaying a well recoverable smooth model (fig. 9)'])

% plot the recovered smooth model
figure(10)
clf
plotconst(malphahat,-pi/2,pi/2);
xlabel('\theta');
ylabel('Intensity');
bookfonts
print -deps c4fsmoothmodre.eps

disp(['Displaying the recovered model from the noisy data associated'...
    ' with the previous plot (fig 10)'])
