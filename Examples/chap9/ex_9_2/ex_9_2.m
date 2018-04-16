% Example 9.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Make x, y, and sigma global.
global x;
global y;
global sigma;

% Generate the data set.
x=(1.0:0.25:7.0)';
ptrue=[1.0; -0.5; 1.0; -0.75];
ytrue=rawfunc(x,ptrue);
sigma=0.01*ones(size(ytrue));
y=ytrue+sigma.*randn(size(ytrue));

dof=length(y)-length(ptrue);

% Output the data set.
[x,y]

% Now, go ahead and solve the problem with many different random starting 
% points.
N=20;
models=zeros(4,N);
chisqvals=zeros(N,1);
normg=zeros(N,1);
itercnts=zeros(N,1);
for i=1:N
  p0=2*(rand(4,1)-0.5*ones(4,1));
  [pest,itercnts(i)]=lm('fun','jac',p0,1.0e-14,500);
  models(:,i)=pest;
  chisqvals(i)=norm((rawfunc(x,pest)-y)./sigma,2)^2;
  ptest(i) = 1-chi2cdf(chisqvals(i),dof);
  normg(i)=norm(jac(pest)'*fun(pest));
end

% With the randomly generated starting points used here, it happens
% that the best of the locally optimal solutions is found from the 
% second starting point, the next best is found from the fifth 
% starting point, and the third best solution is found from
% starting point one.

% Make a table of the locally optimal solutions. (2, 5, 1).  
table=[models(1,2) models(2,2) models(3,2) models(4,2) chisqvals(2) ptest(2);
       models(1,5) models(2,5) models(3,5) models(4,5) chisqvals(5) ptest(5);
       models(1,1) models(2,1) models(3,1) models(4,1) chisqvals(1) ptest(1)]

% Find the best parameters, and covariance matrix.
pbest=models(:,2)
Jbest=jac(pbest);
covbest=inv(Jbest'*Jbest);
corbest=zeros(4,4);
for i=1:4
  for j=1:4
    corbest(i,j)=covbest(i,j)/(sqrt(covbest(i,i)*covbest(j,j))); 
  end 
end
disp('The correlation matrix is');
corbest 

%get the fit data
yfit=rawfunc(x,pbest);

% Plot the fitted model and data points.
figure(1)
clf
plot(x,yfit,'k');
hold on;
% one sigma error bars.
errorbar(x,y,sigma,'ko');
xlabel('x');
ylabel('y');
bookfonts
ylim([0 1.2])

disp('Displaying the fit model and the data (fig. 1)')
print -deps2 c9fnlregfit.eps

% Plot the normalized residuals
figure(2)
clf
plot(x,fun(pbest),'ko');
axis([0 7 -2 2]);
xlabel('x');
ylabel('Normalized Residual');
bookfonts

disp('Displaying the normalized residuals (fig. 2)')
print -deps2 c9fresiduals.eps
