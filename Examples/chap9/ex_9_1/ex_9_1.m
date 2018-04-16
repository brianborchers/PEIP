% Example 9.1
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Global variables
% H     - the recorded head for each time
% TM    - the times the head was recorded
% SIGMA - the standard deviation for a time
% D     - the distance between the wells
% Q     - the volume of the slug
global H;
global TM;
global SIGMA;
global Q;
global D;
%
% Load the data.  Head is measured to the nearest centimeter.  
%
load slugdata.mat
%
% Fixed parameter values.
%
D=60;
Q=50;
% We'll use sigma=1cm.  
SIGMA=0.01*ones(size(H));

% The unknown/estimated parameters are S=p(1) and T=p(2).
p0=[0.001; 1.0];

% Solve the least squares problem with LM.
[pest,itercnt]=lm('fun','jac',p0,1.0e-12,100);

% Check the chi^2 value.
chi2=norm(fun(pest))^2;
pvalue=1-chi2cdf(chi2,length(H)-2);
disp(['The estimate S=',num2str(pest(1)),' T=', num2str(pest(2)), ...
    ' has chi squared ', num2str(chi2), ' which gives a p-value of '...
    num2str(pvalue)])

% Compute the Jacobian, covariance matrix, correlation matrix, and confidence 
% intervals.
Jest=jac(pest);
C=inv(Jest'*Jest);
m=max(size(C));
for i=1:m
  for j=1:m
    Corm(i,j)=C(i,j)/(sqrt(C(i,i))*sqrt(C(j,j)));
  end
end

% display the independent confidence intervals
fprintf(1,'S=%.5f +- %.5f\n',[pest(1) 1.96*sqrt(C(1,1))]);
fprintf(1,'T=%.5f +- %.5f\n',[pest(2) 1.96*sqrt(C(2,2))]);
disp('Covariance matrix is');
C
disp('Correlation matrix is');
Corm

%determine the fitted head
tfit=0.001:0.1:51;
for i=1:length(tfit),
  hfit(i)=Q*exp(-D^2*pest(1)/(4*pest(2)*tfit(i)))/(4*pi*pest(2)*tfit(i));
end

% Now, plot the fitted head
figure(1)
clf
plot(tfit,hfit,'k-');
hold on
errorbar(TM,H,SIGMA,'ko');
xlabel('Time (hours)');
ylabel('Head (m)');
bookfonts

disp('Displaying the predicted head (fig. 1)')
print -deps c9fslugplot.eps

% generate the chi^2 values over a region
[X,Y]=meshgrid(0.0001:0.0001:0.01,0.01:0.01:2.0);
Z=zeros(size(X));
[m,n]=size(X);
for i=1:m
  for j=1:n
    Z(i,j)=norm(fun([X(i,j); Y(i,j)]))^2;
  end
end  

% Produce a contour plot of the chi^2 function.
figure(2);
clf
[CS,HHH]=contour(X,Y,Z,[10 100  1000 3000]);
% make sure that spacing, size and color of the contours are appropriate
clabel(CS,HHH,'LabelSpacing',300,'Fontsize',14);
for i=1:length(HHH)
  set(HHH(i),'EdgeColor','k');
end
hold on
plot([1.83 2.26 2.26 1.83 1.83]*1e-3, [0.54 0.54 0.64 0.64 0.54],...
    'k--','Linewidth',1.5,'Color',[0.5 0.5 0.5]);
ylim([0 2.2])
xlim([0 6.2]*1e-3)
xlabel('Storage Coefficient, S');
ylabel('Trasmissivity, T (m^2/hr)');
bookfonts

disp('Displaying the chi squared contours (fig. 2)')
print -deps c9fslugcontour.eps


% create chi^2 contours nearer the fitted values
[X,Y]=meshgrid(0.0017:0.000025:0.0025,0.50:0.001:0.70);
Z=zeros(size(X));
[m,n]=size(X);
for i=1:m
  for j=1:n
    Z(i,j)=norm(fun([X(i,j); Y(i,j)]))^2;
  end
end  

% compute linearized ellipsoid
[u,lam]=eig(inv(C));
deltachisq=chi2inv(0.95,2);
delta=sqrt(deltachisq);
%generate a vector of angles from 0 to 2*pi
theta=(0:.01:2*pi)';
%calculate the x and y components of the ellipsoid for all angles
r=zeros(length(theta),2);
r(:,1)=(delta/sqrt(lam(1,1)))*u(1,1)*...
    cos(theta)+(delta/sqrt(lam(2,2)))*u(1,2)*sin(theta);
r(:,2)=(delta/sqrt(lam(1,1)))*u(2,1)*...
    cos(theta)+(delta/sqrt(lam(2,2)))*u(2,2)*sin(theta);


% Plot a second contour plot over a parameter range that is much closer to
% the fitted values.

figure(3)
clf
% plot a contour near the preferred solution
[CS,HHH]=contour(X,Y,Z,...
    round(100*[chi2+chi2inv(0.9,2), chi2+chi2inv(0.9,2)])/100);
set(HHH(1),'EdgeColor','k');
clabel(CS,HHH,'LabelSpacing',400,'Fontsize',14);
hold on
% plot a contour farther from the preferred solution
[CS,HHH]=contour(X,Y,Z,...
    round(100*[chi2+chi2inv(0.95,2), chi2+chi2inv(0.95,2)])/100,'Linewidth',3);
set(HHH(1),'EdgeColor','k');
clabel(CS,HHH,'LabelSpacing',200,'Fontsize',14);
% plot a contour even farther from the preferred solution
[CS,HHH]=contour(X,Y,Z,...
    round(100*[chi2+chi2inv(0.99,2), chi2+chi2inv(0.99,2)])/100);
set(HHH(1),'EdgeColor','k');
clabel(CS,HHH,'LabelSpacing',400,'Fontsize',14);

% Put in the individual 95% CI's.
plot([pest(1)-1.96*sqrt(C(1,1)), pest(1)-1.96*sqrt(C(1,1))],[0.5 0.7],'k--');
plot([pest(1)+1.96*sqrt(C(1,1)), pest(1)+1.96*sqrt(C(1,1))],[0.5 0.7],'k--');
plot([0.0017 0.0025],...
    [pest(2)-1.96*sqrt(C(2,2)), pest(2)-1.96*sqrt(C(2,2))],'k--');
plot([0.0017 0.0025],...
    [pest(2)+1.96*sqrt(C(2,2)), pest(2)+1.96*sqrt(C(2,2))],'k--');

% mark the estimated point
plot(pest(1),pest(2),'k.','markersize',14);
text(pest(1),pest(2),' m^{\ast}','Fontsize',20,'Fontweight','bold');

% plot the linearized error ellipsoid
plot(pest(1)+r(:,1),pest(2)+r(:,2),'k--','Color',[0.7 0.7 0.7],'linewidth',3);
hold off

xlabel('Storage Coefficient, S');
ylabel('Transmissivity, T (m^2/hr)');
bookfonts
ylim([0.54 0.64])
xlim([1.83 2.26]*1e-3)

disp(['Displaying contours of the chi squared value, and linearized error'...
    ' ellipse (fig. 3)'])
print -deps c9fslugcontour2.eps
