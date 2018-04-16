% Examples 4.4 and 4.5
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% basic control parameters
noise = 0.0002;
maxdepth=1000;
depth=0:maxdepth;

% generate the true smooth model
strue=slowfunctionsmooth(depth);

% create the noisy observed data
%
% dobs is a vector containing the depths of the different observations.
% dd is difference in depth between successive observations.
% t contains the noisy observations.
%
%forward model borehole arrival time data using function 'slowfunctionsmooth'
i=0;
dd=20;
for d=dd:dd:1000
  i=i+1;
  dobs(i,1)=d;
  t(i,1)=quadl(@slowfunctionsmooth,0,d)+noise*randn();
end


% set up the system matrix
% the number of data points collected
M=length(dobs);
%Number of model parameters
N=M;
%System matrix
G=zeros(M,N);
for i=1:M
  G(i,1:i) = ones(1,i)*dd;
end

figure(1)
clf
plot((0:length(strue)-1),strue*1000,'k');
xlabel('Depth (m)');
ylabel('True Slowness (s/km)');
bookfonts

disp('Displaying the true slowness model (fig. 1)')

%data covariance matrix
%
CD=noise^2*eye(M);
%
% %set up the prior
mprior=linspace(strue(1),strue(end),M)';
%mprior=0.29*1e-3*ones(M,1);
clen=5;
sigkern=2;
%root correlation kernel
ckern=zeros(2*M,1);
%start with a triangle root kernel
ckern(M-clen:M+clen)=bartl(2*clen+1)';
%Gaussian root kernel
%ckern(M-clen:M+clen)=exp((-1/2)*(-clen:clen).^2/sigkern^2);
%cfun, the autocorrelation of ckern, defines the uniform correlation structure
%for the desired matrix
cfun=xcorr(ckern,M,'coeff');
%populate a correlation matrix
CorrM=zeros(M,M);
for i=1:M
    c=circshift(cfun,i+M);
    CorrM(:,i)=c(1:M);
end

%scale the correlation matrix to obtain a uniform covariance matrix for the
%prior
sigmam=2e-5;
CM=CorrM*sigmam^2;

[covmp,mmap]=bayes(G,mprior,CM,t,CD);


figure(2)
plotconstc(mprior*1000,0,maxdepth,'k');
hold on
plot(0:length(strue)-1,strue*1000,'k','Color',[.7 .7 .7],'LineWidth',3);
plotconstc(mprior*1000-1.96*1000*sqrt(diag(CM)),0,maxdepth,'k-.');
plotconstc(mprior*1000+1.96*1000*sqrt(diag(CM)),0,maxdepth,'k-.');
hold off
xlabel('Depth (m)');
ylabel('Slowness (s/km)');
bookfonts
print -depsc2 c11fvspprior.eps

disp('Displaying the prior slowness distribution (fig. 2)')

%MMAP solution
figure(3)
clf
plotconst(mmap*1000,0,maxdepth);
hold on
plot(0:length(strue)-1,strue*1000,'k','Color',[.7 .7 .7],'LineWidth',3);
plotconstc(mmap*1000-1.96*1000*sqrt(diag(covmp)),0,maxdepth,'k--');
plotconstc(mmap*1000+1.96*1000*sqrt(diag(covmp)),0,maxdepth,'k--');
hold off
ylim([0.2 0.36])
xlabel('Depth (m)');
ylabel('Slowness (s/km)');
bookfonts

ylabel('Slowness (s/km)');
xlabel('Depth (m)');
ylim([0.2 0.36])
print -depsc2 c11fvsp1bayes.eps
disp('Displaying the first MAP model (fig. 3)')

figure(4)
plotconst(CorrM(25,:),-500,500)
xlim([-500 500])
ylabel('a_i')
xlabel('Lag (m)')
bookfonts
print -deps2 c11fvsp1corr.eps

disp('Displaying the correlation function for the first MAP model (fig. 4)')

%set up a smoother prior
mprior=linspace(strue(1),strue(end),M)';
clen=10;
sigkern=2;
%root correlation kernel
ckern=zeros(2*M,1);
%start with a triangle root kernel
ckern(M-clen:M+clen)=bartl(2*clen+1)';
%Gaussian root kernel
%ckern(M-clen:M+clen)=exp((-1/2)*(-clen:clen).^2/sigkern^2);
%cfun, the autocorrelation of ckern, defines the uniform correlation structure
%for the desired matrix
cfun=xcorr(ckern,M,'coeff');
%populate a correlation matrix
CorrM=zeros(M,M);
for i=1:M
    c=circshift(cfun,i+M);
    CorrM(:,i)=c(1:M);
end

%scale the correlation matrix to obtain a uniform covariance matrix for the
%prior
sigmam=2e-5;
CM=CorrM*sigmam^2;

[covmp,mmap]=bayes(G,mprior,CM,t,CD);

%MMAP solution
figure(5)
clf
plotconst(mmap*1000,0,maxdepth);
hold on
plot(0:length(strue)-1,strue*1000,'k','Color',[.7 .7 .7],'LineWidth',3);
plotconstc(mmap*1000-1.96*1000*sqrt(diag(covmp)),0,maxdepth,'k--');
plotconstc(mmap*1000+1.96*1000*sqrt(diag(covmp)),0,maxdepth,'k--');
hold off
ylabel('Slowness (s/km)');
xlabel('Depth (m)');
bookfonts;
print -depsc2 c11fvsp2bayes.eps

disp('Displaying the second MAP model (fig. 5)')

figure(6)
plotconst(CorrM(25,:),-500,500)
xlim([-500 500])
xlabel('Lag (m)')
ylabel('a_i')
bookfonts
print -deps2 c11fvsp2corr.eps

disp('Displaying the correlation function for the second MAP model (fig. 6)')
