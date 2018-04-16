% Example 10.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% Demonstrate the adjoint method for computing the derivative of a solution
% of the time dependent heat equation
%
%   partial(u)/partial(t)=D partial^2(u)/partial^2(x)
%   u(0,t)=0, u(1,t)=0
%
% Given u(xi,T) at a set of points x1, x2, ..., xk, we'd like to
% estimate u(x,0).  
%
clear
%
% Reset the RNG.
%
rand('seed',0);
randn('seed',0);
%
% Global variables.
%
global m;
global n;
global deltax;
global x;
global deltat;
global xpoints;
global d;
global D;
global Fvec;
global alpha;
global L;
%
% Setup model discretization
%
n=1000;
xmin=-1;
xmax=1;
deltax=0.002;
x=(deltax:deltax:(n*deltax))'-1;

deltat=1.0e-3;
deltat=deltat*1e4*5;
%
% Set the thermal diffusion coefficient (m^2/s).
%
D=1e-6;
%
% Set the time period T s.
%
T=1e4;
m=T/deltat;
%
% Initial condition (true solution, 0<=x<=1.
%
u0true=zeros(size(x));
u0true(500:600)=10;
u0true(200:250)=5;
u0true(775:800)=3;
%
% Compute the solution.
%
[uT,A,B]=forward(m,n,deltat,deltax,D,u0true);
%
% Get the data points.
%
xpoints=(25:25:975)';
d=uT(xpoints);
%
% Add noise to the data.
%
noiselevel=0.1;
d=d+noiselevel*randn(size(d));
%
% Plot the solution at time T and the noisy data points. 
%
figure(1);
clf
plot(x,uT,'k');
hold on
plot(x(xpoints),d,'ko');
xlabel('x (m)')
ylabel('\Delta T (^oK)');
bookfonts
print -deps2 c10fheatdata.eps

%
% Now, setup and solve the inverse problem. 
%
%
% First, setup the regularization
%
%zeroth-order Tikhonov matrix for the objective function
L=speye(n);

%search over a range of regularization parameters
alphas = logspace(-5,-2.5,20);
for i=1:length(alphas)
    alpha=alphas(i);
%
% Setup a first guess for u0(a uniform temperature perturbation)
%
u0guess=10*ones(n,1);
%
% Call the conjugate gradient method to find the best solution.
%
[u0min,fu0min]=conjg(@objfun,@grad,u0guess,1.0e-15);
%
% Plot the solution and compare with u0true.
%
figure(3);
clf
plot(x,u0true,'k');
hold on
plot(x,u0min,'k--');
legend('u0 true',['recovered model for \alpha = ',num2str(alphas(i))]);
xlabel('x (m)')
ylabel('\Delta T (^oK)');
bookfonts;

% Plot u(x,T) corresponding to that u0min.
%
figure(4);
clf
uTmin=forward(m,n,deltat,deltax,D,u0min);
plot(x,uTmin,'k');
hold on
plot(x(xpoints),d,'ko');
xlabel('x (m)')
ylabel('\Delta T (^oK)');
bookfonts

%discrepancy criterion evaluation
disp(['Norm of Residual vector for alpha =',num2str(alphas(i))])
resid(i)=norm(d-uTmin(xpoints));
u0save(:,i)=u0min;
mnorm(i)=norm(u0min);

end

discrep=sqrt(length(xpoints))*noiselevel;
disp(['Discrepancy Principle Residual Criterion: ',num2str(discrep)])

figure(5)
clf
plot(resid,mnorm)
hold on
plot(resid,mnorm,'o','markersize',10)
for i=3:2:length(alphas)
    text(resid(i)+0.03,mnorm(i),num2str(alphas(i),'%6.5f'),'fontsize',20);
end
plot([discrep,discrep],[min(mnorm),max(mnorm)],'--','linewidth',3)
text(discrep+0.03,min(mnorm)+5,['\delta = ',num2str(discrep)],'fontsize',25)
xlabel('Residual Norm || G(m) - d ||_2')
ylabel('Solution Norm ||m||_2')
bookfonts
hold off
axis tight
print -deps2 c10fheattradeoff.eps

%with above search, the closest model to meeting the discrepancy principle
%is model number 11
figure(6);
clf
plot(x,u0true,'k');
hold on
plot(x,u0save(:,11),'k--');
legend('u0 true',['recovered model for \alpha = ',num2str(alphas(11))]);
xlabel('x (m)')
ylabel('\Delta T (^oK)');
bookfonts
ylim([-2 13])
print -deps2 c10fheatsoln.eps

