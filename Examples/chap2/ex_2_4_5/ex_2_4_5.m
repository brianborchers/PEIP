% Examples 2.4 and 2.5 
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

%make sure we have a clean environment
clear;
rand('state',0);
randn('state',0);

%
% The code that generated the data.
%
%beta0 = 10;
%beta1 = 100;
%beta2 = 9.8;
%for t=1:10
%x(t,1)=t;
%sigma(t,1)=8;
%y(t,1)=beta0+beta1*t-(1/2*beta2)*t^2+sigma(t,1)*randn;
%end

% Load precomputed data
load data1.mat
t=data1(:,1);
y=data1(:,2);
sigma=data1(:,3);
N=length(t);

% Introduce the outlier
y(4)=y(4)-200;

disp('displaying t,y,sigma')
[t, y, sigma ]

% Build the parabolic system matrix

G = [ ones(N,1), t, -1/2*t.^2 ];

% Apply the weighting

yw = y./sigma;
Gw = G./[sigma,sigma,sigma];

% Solve for the least-squares solution
disp('least-squares solution')
m2 = Gw\yw


% Solve for the 1-norm solution
disp('1-norm solution')
m1 = irls(Gw,yw,1.0e-5,1.0e-5,1,125)


% Get the covariance matrix
ginv = inv(Gw'*Gw)*Gw';

disp('Covariance matrix')
covm = ginv*ginv'


% Get the 1.96-sigma (95%) conf intervals
disp('95% parameter confidence intervals (m-, mest, m+) on 2-norm solution')
del = 1.96*sqrt(diag(covm));
[m2-del , m2 , m2+del]


%
% Because there are 3 parameters to estimate, we have N-3 degrees
% of freedom.
%
dof = N-3;
disp(['Chi-square misfit for ',num2str(dof),' dof'])
chi2 = norm((y - G*m2)./sigma)^2

% Find the p-value for this data set
disp('Chi-square p-value')
p = 1-chi2cdf(chi2,dof)


% Find the parameter correlations
s=sqrt(diag(covm))
disp('Parameter correlations for 2-norm solution')
r = covm./(s*s')


% Compute the predicted data from the two models
tt=min(t)-1:.05:max(t)+1;
mm1=m1(1)+tt*m1(2)-.5*tt.^2*m1(3);
mm2=m2(1)+tt*m2(2)-.5*tt.^2*m2(3);

% Plot the models and their data fit to the data
figure(1)
clf
plot(tt,mm1,'k')
hold on
plot(tt,mm2,'--k');
errorbar(t,y,sigma,'ok');
xlabel('Time (s)');
ylabel('Elevation (m)');       
bookfonts
legend('L_1 Fit','L_2 Fit','Location','SouthEast');
hold off

disp('displaying data vs. model fits (fig 1)');
print -deps2 c2fl1example.eps

% Monte Carlo Section for the 1-norm model
% The 1-norm estimated data
y0 = G*m1; 

nreal=1000;

disp('generating Monte Carlo realizations');
%generate the nreal Monte Carlo solutions
for j = 1:nreal
  % Generate a trial data set of perturbed, weighted data
  ytrial = y0+sigma.*randn(N,1);
  ywtrial=ytrial./sigma;
  
  % Store the 1-norm parameters and misfit
  mmc(j,:)=irls(Gw,ywtrial,1.0e-5,1.0e-5,1,500)';
  chimc(j)= norm((G*mmc(j,:)'-y0)./sigma)^2;
end

% figure out how many Monte Carlo points are in the %95 confidence region
disp('confidence interval inclusion check (should be about 0.95):') 
nnz(chimc<=chi2inv(.95, N-3))/nreal

% Calculate the covariance of the parameters in the realizations
disp('Emperical covariance of m1 models');
covmemp=mmc-ones(nreal,1)*mean(mmc);
covmemp=(covmemp'*covmemp)/nreal

% Get the 1.96-sigma (95%) conf intervals
disp('95% parameter confidence intervals (m-, mest, m+) on 1-norm solution')
del = 1.96*sqrt(diag(covmemp));
[m1-del , m1 , m1+del]

% Plot a histogram of the chi-square values
figure(2)
clf
hist(chimc);
title('1000 Monte Carlo Chi-square Values')
ylabel('N')
xlabel('\chi^2')
bookfonts
disp('displaying 1000 Monte Carlo chi-square values (fig 2)');

% Plot the histogram for each model element for all realizations
figure(3)
clf

subplot(1,3,1)
hist(mmc(:,1))
title('m_1 (m)')
bookfonts

subplot(1,3,2)
hist(mmc(:,2))
title('m_2 (m/s)')
bookfonts

subplot(1,3,3)
hist(mmc(:,3))
title('m_3 (m/s^2)')
bookfonts

disp('displaying histogram of 1000 monte-carlo models (fig 3)');

% Calculate the 95% error region for each pair of parameters 
% first get the correct covariance matrix for each pair
vp1=covmemp(1:2,1:2);
vp2=covmemp([1 3],[1 3]);
vp3=covmemp(2:3,2:3);

% get axis directions and lengths for each
[u1,lam1]=eig(inv(vp1));
[u2,lam2]=eig(inv(vp2));
[u3,lam3]=eig(inv(vp3));

% generate the data to plot
% Note that we use 2 degrees of freedom here, because we're
% constructing confidence ellipsoids for 2 parameters at a time.
%
theta = 0:.01:2*pi;
delta = sqrt(chi2inv(.95, 2));
% the radii for each ellipsoid
r1=delta./sqrt(diag(lam1));
r2=delta./sqrt(diag(lam2));
r3=delta./sqrt(diag(lam3));

%compute the ellipsoids
e1=r1(1)*u1(:,1)*cos(theta)+r1(2)*u1(:,2)*sin(theta);
e2=r2(1)*u2(:,1)*cos(theta)+r2(2)*u2(:,2)*sin(theta);
e3=r3(1)*u3(:,1)*cos(theta)+r3(2)*u3(:,2)*sin(theta);

%move then to be centered on the estimated parameters
e1=e1+[m1(1); m1(2)]*ones(1,length(theta));
e2=e2+[m1(1); m1(3)]*ones(1,length(theta));
e3=e3+[m1(2); m1(3)]*ones(1,length(theta));

% Plot the model elements and 95% confidence region for each pair of elements
figure(4)
clf

% The m1, m2 ellipsoid
subplot(1,3,1)
plot(mmc(:,1),mmc(:,2),'b*')
hold on
plot(e1(1,:),e1(2,:),'k')
xlabel('m_1 (m)')
ylabel('m_2 (m/s)')
bookfonts

% The m1, m2 ellipsoid
subplot(1,3,2)
plot(mmc(:,1),mmc(:,3),'b*')
hold on
plot(e2(1,:),e2(2,:),'k')
xlabel('m_1 (m)')
ylabel('m_3 (m/s^2)')
bookfonts

% The m2, m3 ellipsoid
subplot(1,3,3)
plot(mmc(:,2),mmc(:,3),'b*')
hold on
plot(e3(1,:),e3(2,:),'k')
xlabel('m_2 (m/s)')
ylabel('m_3 (m/s^2)')
bookfonts
disp('displaying projections of monte-carlo models (fig 4)');

%q-q plots of parameters, demonstrating their normal distributions
figure(5)
clf


subplot(3,1,1)
qqplot(mmc(:,1))
ylabel('Quantiles, m_1')
xlabel('Quantiles, Standard Normal')
bookfonts


subplot(3,1,2)
qqplot(mmc(:,2))
ylabel('Quantiles, m_2')
xlabel('Quantiles, Standard Normal')
bookfonts


subplot(3,1,3)
qqplot(mmc(:,3))
ylabel('Quantiles, m_3')
xlabel('Quantiles, Standard Normal')
bookfonts
disp('displaying QQ plots of monte-carlo models (fig 5)');
