% Example 6.3 part b
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% This script produces the L-curve and GCV function for Example 6.3 for

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% set up variables we will want to be global
global d;
global GT;
global LT;
global G;
global L;
global alpha;

% Load in the problem data.  
load ex_6_3_30000.mat

% Setup G' and L';

% Load in the problem data.  
load ex_6_3_30000.mat

% Setup G' and L';
GT=G';
LT=L';

% This took about 5 hours on a Phenom II 2.8 GHz processor, for example
% Setup a range of values of the regularization parameter alpha.
alphas=[0.01; 0.03; 0.05; 0.07; 0.1; 0.3; 0.5; 0.7; 1; 3;  5; ...
	7; 10; 30; 50; 70;  100];

% Allocate space for the function outputs.
normgmd=zeros(size(alphas));
normlm=zeros(size(alphas));
v0=zeros(size(alphas));

% Loop through the values of alpha, computing the GCV function for each one.
for k=1:length(alphas);
  % Compute the GCV function (v0), norm(G*mreg-d), and norm(L*mreg).
  [v0(k),num,denom,normgmd(k),normlm(k)]=gcv(alphas(k));

  % Print out the results for this point.
  fprintf('alpha=%f, norm(G*m-d)=%f, norm(L*m)=%f, G(alpha)=%f\n',...
	  [alphas(k); normgmd(k); normlm(k); v0(k)]);
end

% Plot the L-curve.
figure(1)
clf
plot(normgmd,normlm,'k-');
xlabel('Residual Norm ||Gm - d||_{2}');
ylabel('Solution Seminorm ||Lm||_{2}');
bookfonts
xlim([0 80])
hold on
plot(normgmd(7),normlm(7),'ko','markersize',14);
hold off
print -deps2 c6flstoch.eps

% Plot the GCV curve
figure(2)
clf
loglog(alphas,v0,'k-');
hold on
loglog(alphas(7),v0(7),'ko','markersize',14);
hold off
xlabel('Regularization Parameter \alpha');
ylabel('g(\alpha)');
bookfonts
print -deps2 c6fgcvstoch.eps
% save what will be needed for the next example script (ex_6_3c)
save workspace_ex_6_3b.mat
