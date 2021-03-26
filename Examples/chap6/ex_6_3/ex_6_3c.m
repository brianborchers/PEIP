% Example 6.3 part b
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% Resolution portion.  This script computes the regularized solution for 
% alpha=0.5, estimates the diagonal of the resolution matrix, and then 
% compares the estimate against accurately computed values in 100 randomly
% chosen columns.

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% set up the variables we will want to be global
global GT;
global G;
global L;
global LT;
global d;
global alpha;

% Load in the problem data.
load ex_6_3_30000.mat

% Calculate some stuff not stored in the file.
GT=G';
LT=L';

% Some size parameters.
[m,n]=size(G);
[r,s]=size(L);

% Pad d with extra 0's.
dpad=sparse([d; zeros(r,1)]);

% Set the regularization parameter.
alpha=0.5;

disp('This script took about 1.25 hours on a Phenom II 2.8 GHz processor');
disp('finding the regularized solution.');
% Find the regularized solution.
[mreg,FLAG]=lsqr(@Gregmult,dpad,[],1500);

% Output some measures of the result.
disp(['norm(d) is ', num2str(norm(d))]);
disp(['norm(G*m-d) is ', num2str(norm(G*mreg-d))]);
disp(['norm(alpha*L*m) is ', num2str(norm(alpha*(L*mreg)))]);

% Now, repeat the process of estimating the diagonal ntimes.
ntimes=10;
diags=zeros(ntimes,n);
for k=1:ntimes
  % In this case, we'll just use 256 random vectors.
  numcols=256;
  fprintf('Starting diagonal estimation %d\n',k);

  % Estimate the diagonal of R.
  diagR=diagestrand(@Rmult,n,numcols);

  % Store the estimate.  
  diags(k,:)=diagR';
end

% Find the median of the estimates
diagR=median(diags);

% Now, check a random collection diagonal elements to see how well we did. 
% This is done by picking random nonzero columns and computing the resolution 
% of a spike at the corresponding model element. The element corresponding to 
% the random column is what would on the diagonal.
% Note that a zero column has no rays through it, so it makes sense that these 
% can never be resolved.
ntests=100;
cols=zeros(ntests,1);
truevals=zeros(ntests,1);
ests=zeros(ntests,1);
relerr=zeros(ntests,1);
abserr=zeros(ntests,1);
radius=zeros(ntests,1);
raydensity=zeros(ntests,1);

disp('Calculating 100 real diagonal values of the resolution matrix');
for i=1:ntests
  % Pick a random column of G/random element from m
  l=floor(rand*n)+1;
  % Skip any columns that are 0, since the corresponding column of R is also 0.
  while (norm(G(:,l))==0)
    l=floor(rand*n)+1;
  end

  % Finally got a non zero column.  
  % Now, compute this column of R by inverting the spike.
  m=zeros(n,1);
  m(l)=1;
  rhs=G*m;
  rhspad=sparse([rhs; zeros(r,1)]);
  [ker,FLAG]=lsqr(@Gregmult,rhspad,[],5000);

  % Now, compare the actual value of R(l,l) to the estimate.
  cols(i)=l;
  truevals(i)=ker(l);
  ests(i)=diagR(l);
  relerr(i)=abs(ker(l)-diagR(l))/abs(ker(l));
  abserr(i)=abs(ker(l)-diagR(l));

  %find the radius of this model element center from the center of the
  %model
  [msize,~,~]=size(mtrue);
  [I,J,K]=ind2sub(size(mtrue),l);
  radius(i)=sqrt((I-msize/2)^2+(J-msize/2)^2+(K-msize/2)^2);
  fprintf(1,['estimated r(%d)=%.3f, accurate r(%d)=%.3f, relerr=%.2f,'...
      ' abserr=%.3f element radius=%.2f\n'], ...
	    [l; diagR(l); l; ker(l); relerr(i); abserr(i); radius(i)]);

  % While we're at it, compute the ray density for this column.
  raydensity(i)=nnz(G(:,l));
end


% Print out the mean absolute error over the 100 randomly chosen columns.
disp(['The mean absolute error was ', num2str(mean(abserr))]);

% Plot a scatter plot showing the estimates of the resolutions
% versus the accurately computed values.

figure(1)
clf
plot(truevals,ests,'ko');
xlabel('Computed Values of R_{m} Diagonal');
ylabel('Stochastic Estimates of R_{m} Diagonal');
bookfonts
print -deps2 c6fstochasticest.eps

% Plot a scatter plot showing the ray density vs. the accurately
% computed resolution values.

figure(2)
clf
bookfonts
plot(truevals,raydensity,'ko');
xlabel('Computed Values R_{m} Diagonal');
ylabel('Ray Count');
print -deps2 c6fraydensity.eps

% Plot a scatter plot showing the ray density vs. the distance from the middle 
% of the region

figure(3)
clf
plot(raydensity,radius,'ko');
xlabel('Ray densities');
ylabel('Radius of Element From Center of Model');
bookfonts

% plot a scatter plot showing the accurate resolution diagonals vs. the 
% distance from the middle of the region

figure(4)
clf
plot(radius,truevals,'ko');
xlabel('Radius of Element From Center of Model');
ylabel('Computed Values of R_{m} Diagonal');
bookfonts

% save the results for later analysis if desired
%save workspace_ex_6_3c.mat

