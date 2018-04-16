% Example 2.3
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Generate the x and y values
x=(25:3:97)';
ytrue=10*x;
N=length(x);

% compute the noisy y values
y=ytrue+0.05*randn(size(x)).*ytrue;

% Do the regression.  
G=[ones(size(x)) x];
m=inv(G'*G)*G'*y;

% display the model.
disp('the model parameters are');
m

% Compute the modeled y's.
ymod=G*m;

% Compute the residuals and estimate s and C.
r=ymod-y;
s=sqrt(norm(r)^2/(N-2))
C=s^2*inv(G'*G)

% Compute confidence intervals for the parameters.
disp('half widths')
hw=tinv(0.975,N-2)*sqrt(diag(C))

% Plots
%  1. Data and fitted line.
%  2. residuals.
figure(1)
clf
plot(x,y,'ko');
hold on
plot(x,ymod,'k');
xlabel('x');
ylabel('y');
bookfonts

disp('Displaying Data and Linear Regression Line (fig 1)')
print -deps2 c2ffitted.eps

figure(2)
clf
plot(x,y-ymod,'ko');
xlabel('x');
ylabel('r');
bookfonts

disp('Displaying Model Residuals vs. x (fig 2)');
print -deps2 c2fresid.eps

% Now, look at the scaled problem.
% W has the inverse of the y values on the diagonal
W=inv(diag(y));
Gw=W*G;
yw=W*y;
mw=inv(Gw'*Gw)*Gw'*yw;

% display the model.
disp('the model parameters are');
mw

% Compute the modeled weighted y's.
ymodW=Gw*mw;

% Compute the residuals and estimate s and C.
rw=ymodW-yw;
sw=sqrt(norm(rw)^2/(N-2))
Cw=sw^2*inv(Gw'*Gw)

% Compute confidence intervals for the parameters.
disp('half widths')
hww=tinv(0.975,N-2)*sqrt(diag(Cw))

% Compute the modeled unweighted y's.
ymodw=G*mw;

% Plots
%  3. Data and fitted line.
%  4. residuals.
figure(3)
clf
plot(x,y,'ko');
hold on
plot(x,ymodw,'k');
xlabel('x');
ylabel('y');
bookfonts

disp('Displaying Data and Weighted Linear Regression Line (fig 3)');
print -deps2 c2ffittedsc.eps

figure(4)
clf
plot(x,rw,'ko');
xlabel('x');
ylabel('r');
bookfonts

disp('Displaying Weighted Model Residuals vs. x (fig 4)');
print -deps2 c2fresidsc.eps
