% Example 4.7
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% (originally suggested by Gary Pavlis, Indiana University).

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

%load particular realization of the data vector used in book
load vsp

%Number of model parameters
k=1;
N=k*M;
dd=20/k;

%System matrix
G=zeros(M,N);
for i=1:M
  G(i,1:k*i) = ones(1,k*i)*dd;
end

% Get least squares solution
ngcvpoints=1000;

% Get roughening matrix and generalized svd (using cgsvd.m)
L1 = get_l_rough(N,1);
[U1,V1,X1,Lam1,MU1]=gsvd(G,L1);
lam=sqrt(diag(Lam1'*Lam1));
mu=sqrt(diag(MU1'*MU1));
p=rank(L1);
sm1=[lam(1:p),mu(1:p)];

% get the gcv values varying alpha
[alpha1,g1,reg_param1]=gcval(U1,sm1,t,ngcvpoints);

% get the minimum value of of the g functions
[ming1,~] = min(g1);

% Plot first order GCV function.
figure(1)
clf
loglog(reg_param1,g1,'-k');
xlabel('\alpha');
ylabel('g(\alpha)');
bookfonts
ax=axis;
hold on
H=loglog(alpha1,ming1,'ok',[alpha1,alpha1],[ming1/ngcvpoints,ming1],'-k');
set(H,'markersize',16);
hold off
axis(ax);
xlim([1 1e4]);
ylim([5e-10 2e-8]);
print -deps2 c4fgcv1curve.eps

disp('Displaying GCV L-curve for first order regularization (fig. 1)')

% Second order Tikhonov solution using GCV
% Get roughening matrix and generalized svd
L2 = get_l_rough(N,2);
[U2,V2,X2,Lam2,MU2]=gsvd(G,L2);
lam=sqrt(diag(Lam2'*Lam2));
mu=sqrt(diag(MU2'*MU2));
p=rank(L2);
sm2=[lam(1:p),mu(1:p)];

% get the gcv values varying alpha
[alpha2,g2,reg_param2]=gcval(U2,sm2,t,ngcvpoints);

%get the minimum value of of the g functions
[ming2,~] = min(g2);

% Plot second order GCV function. 
figure(2)
clf
loglog(reg_param2,g2,'-k');
xlabel('\alpha');
ylabel('g(\alpha)');
bookfonts
ax=axis;
hold on
H=loglog(alpha2,ming2,'ok',[alpha2,alpha2],[ming2/ngcvpoints,ming2],'-k');
set(H,'markersize',16);
hold off
axis(ax);
xlim([1 1e4]);
ylim([7e-10 2e-9]);
print -deps2 c4fgcv2curve.eps

disp('Displaying GCV L-curve for second order regularization (fig. 2)')

% compute the model results
m1=(G'*G+alpha1^2*(L1'*L1))\G'*t;
m2=(G'*G+alpha2^2*(L2'*L2))\G'*t;

% Plot first order Tikhonov regularized GCV solution
figure(3)
clf
plotconst(m1*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'k--');
hold off
xlabel('Depth (m)')
ylabel('Slowness (s/km)')
bookfonts

disp('Displaying GCV solution with first order regularization (fig. 3)')
print -deps2 c4fgcv1model.eps

% Plot second order Tikhonov regularized GCV solution
figure(4)
clf
plotconst(m2*1000,0,maxdepth);
hold on
plot(depth,strue*1000,'k--');
hold off
xlabel('Depth (m)')
ylabel('Slowness (s/km)')
bookfonts
print -deps2 c4fgcv2model.eps

disp('Displaying GCV solution with second order regularization (fig. 4)')

% display the best alphas
disp(['1st order alpha is:  ',num2str(alpha1)]);
disp(['2nd order alpha is:  ',num2str(alpha2)]);

% compute and print the misfits
dinterp=10:20:990;
strueresamp=spline(depth,strue,dinterp)';

disp(['1st order model 2-norm misfit (s/km): ',num2str(1000*norm(m1-strueresamp))]);
disp(['2nd order model 2-norm misfit (s/km): ',num2str(1000*norm(m2-strueresamp))]);
