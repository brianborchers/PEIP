% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

%Global variables for use by the mcmc function calls
global x;
global y;
global sigma;
global stepsize;

disp('Note: This example will take several minutes to run on most computers')

% Generate the data set.
x=(1:0.25:7.0)';
mtrue=[1.0; -0.5; 1.0; -0.75];
ytrue=fun(mtrue,x);
sigma=0.01*ones(size(ytrue));
y=ytrue+sigma.*randn(size(ytrue));

%set the MCMC parameters
%number of skips to reduce autocorrelation of models
skip=10000;
%burn-in steps
BURNIN=100000;
%number of posterior distribution samples
N=4100000;
%MVN step size
stepsize = 0.005*ones(4,1);

% We assume flat priors here
%% initialize model at a random point on [-1,1]
m0=(rand(4,1)-0.5)*2;

[mout,mMAP,pacc]=mcmc('logprior','loglikelihood','generate','logproposal',m0,N);
 
%
% plot autocorrelations.
%

figure(1)
clf
%plot parameter correlations
laglen=10000;
lags=(-laglen:laglen)';
acorr=zeros(2*laglen+1,4);

disp('Calculating autocorrelations for fig. 1')
for i=1:4
      acorr(:,i)=calc_corr(mout(i,:)',laglen);
end

for i=1:4
  subplot(4,1,i);
  bookfonts;
  plot([0 laglen],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',3);
  hold on
  plot(lags(laglen+1:200:laglen*2+1),acorr(laglen+1:200:laglen*2+1,i),'ko');
  hold off
  ylabel(['A ( m_',num2str(i),')'])
  ylim([-0.5 1])
  if i~=4
    set(gca,'Xticklabel',[]);
  end
  bookfonts
end
xlabel('Lag')
bookfonts

print -depsc2 c11MCMCmcorrbefore.eps
disp(['Displaying the autocorrelation of individual parameters before' ...
      ' thinning (fig. 1)']);

%  %track accepted new models to report acceptance rate
%  nacc=0;
%   
%   %sample the posterior here
%   for t = 2:N,
%       %calculate a candidate model
%    c=m(:,t-1)+randn(4,1).*stepsize;
%       %calculate acceptance probability for c,m(:,t-1)
%       lnprobacc=getlnprobacc(m(:,t-1),c,x,y,sigma);
%       if log(rand) < lnprobacc
%           m(:,t)=c;
%           nacc=nacc+1;
%       else
%           m(:,t)=m(:,t-1);
%       end
%   end
disp(['Acceptance Rate: ',num2str(pacc)]);
  
%downsample results to reduce correlation
k=(BURNIN:skip:N);
 
mskip=mout(:,k);
  
%histogram results, and find the modes of the subdistributions as an
%estimte of the MAP model
disp(['m_map','  m_true'])
[mMAP,mtrue]
 
% estimate the 95% credible intervals
for i=1:4
  msort=sort(mskip(i,:));
  m2_5(i) = msort(round(2.5/100*length(mskip)));
  m97_5(i) =  msort(round(97.5/100*length(mskip)));
  disp(['95% confidence interval for m', num2str(i),' is [', num2str(m2_5(i)),',', num2str(m97_5(i)),']'])
end

%plot a scatter plot and histogram of the posterior distribution
mlims=[0.7 1.2; -1 -0.4 ; 0.75 1.5 ; -0.9 -0.6];

figure(2)
clf
for i=1:4
  for j=1:4
    subplot(4,4,4*(i-1)+j)
    if i==j
      hist(mskip(i,:));
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','k')
      set(gca,'Yticklabel',[]);
      xlim(mlims(i,:));
      bookfonts
    else
      plot(mskip(j,:),mskip(i,:),'k.','Markersize',6,'Color',[0.6 0.6 0.6]);
      hold on
      % plot the true answer as a large black dot
      plot(mtrue(j),mtrue(i),'k.','Markersize',24);
      % plot the accepted answers as gray dots
      plot(mMAP(j),mMAP(i),'ko','Markersize',12,'LineWidth',3);
      % plot the 95% ci as a box
      plot([m2_5(j),m97_5(j)],[m2_5(i),m2_5(i)],'k-','LineWidth',1);
      plot([m2_5(j),m97_5(j)],[m97_5(i),m97_5(i)],'k-','LineWidth',1);
      plot([m2_5(j),m2_5(j)],[m2_5(i),m97_5(i)],'k-','LineWidth',1);
      plot([m97_5(j),m97_5(j)],[m2_5(i),m97_5(i)],'k-','LineWidth',1);
      xlim(mlims(j,:));
      ylim(mlims(i,:));
      bookfonts
      hold off
    end
  end
end

subplot(4,4,1)
ylabel('m_1')
bookfonts
subplot(4,4,5)
ylabel('m_2')
bookfonts
subplot(4,4,9)
ylabel('m_3')
bookfonts
subplot(4,4,13)
ylabel('m_4')
bookfonts
subplot(4,4,13)
xlabel('m_1')
bookfonts
subplot(4,4,14)
xlabel('m_2')
bookfonts
subplot(4,4,15)
xlabel('m_3')
bookfonts
subplot(4,4,16)
xlabel('m_4')
bookfonts

print -depsc2 c11MCMCscat.eps
disp('Displaying the true parameters, thinned points, and 95% confidence intervals (fig. 2)');

%plot parameter sample histories
figure(3)
clf
for i=1:4
  subplot(4,1,i)
  plot([1 length(mskip)],[mtrue(i) mtrue(i)],'Color',[0.6 0.6 0.6],'LineWidth',3);
  hold on
  plot(mskip(i,:),'ko')
  hold off
  if i~=4
    set(gca,'Xticklabel',[]);
  end
  xlim([1 length(mskip)])
end
xlabel('Sample Number')
subplot(4,1,1)
ylabel('m_1')
bookfonts
subplot(4,1,2)
ylabel('m_2')
bookfonts
subplot(4,1,3)
ylabel('m_3')
bookfonts
subplot(4,1,4)
ylabel('m_4')
bookfonts

print -depsc2 c11MCMCmhist.eps
disp('Displaying the true parameters and thinned sample histories (fig. 3)');

%plot parameter correlations
figure(4)
clf
laglen=50;
lags=(-laglen:laglen)';
acorr=zeros(2*laglen+1,4);
for i=1:4
  acorr(:,i)=calc_corr(mskip(i,:)',laglen);
  subplot(4,1,i);
  plot([0 laglen],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',3);
  hold on
  plot(lags(laglen+1:laglen*2+1),acorr(laglen+1:laglen*2+1,i),'ko');
  hold off
  ylabel(['A ( m_',num2str(i),')'])
  ylim([-0.5 1])
   bookfonts
  if i~=4
    set(gca,'Xticklabel',[]);
  end
end
xlabel('Lag')
bookfonts
print -depsc2 c11MCMCmcorr.eps
disp('Displaying the autocorrelation of thinned parameter samples (fig. 4)');