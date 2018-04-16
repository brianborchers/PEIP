% Example 6.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Load the raw image in.
img=double(imread('image.png'));

% Build the G matrix.
G=blur(200,3,15);

% Compute the blurred image.
d=G*reshape(img,40000,1);

% Add noise as appropriate.
dn=d+2.0e-2*randn(size(d));

% Use CGLS (with no explicit regularization) to deblur.
[X,rho,eta]=cgls(G,dn,200);


% Plots with no explicit regulariztion.

% Plot raw image
figure(1)
clf
imagesc(img);
colormap(gray);
H=gca;
set(H,'XTick',[]);
set(H,'YTick',[]);
bookfonts

disp('Displaying raw image (fig. 1)')

% Plot blurred image with noise
figure(2)
clf
imagesc(reshape(dn,200,200));
colormap(gray);
H=gca;
set(H,'XTick',[]);
set(H,'YTick',[]);

disp('Displaying blurred image with noise (fig. 2)')
bookfonts
print -deps c6fblur.eps

% Plot CGLS solution with 30 iterations 
figure(3)
clf
imagesc(reshape(X(:,30),200,200));
colormap(gray);
H=gca;
set(H,'XTick',[]);
set(H,'YTick',[]);
bookfonts

disp('Displaying CGLS solution after 30 iterations (fig. 3)')
print -deps c6fblur30.eps

% Plot CGLS solution with 100 iterations 
figure(4)
clf
imagesc(reshape(X(:,100),200,200));
colormap(gray);
H=gca;
set(H,'XTick',[]);
set(H,'YTick',[]);
bookfonts

disp('Displaying CGLS solution after 100 iterations (fig. 4)')
print -deps c6fblur100.eps

% Plot of model norm vs. residual norm for CGLS solution
figure(5)
clf
bookfonts
loglog(rho(1:200),eta(1:200),'ko');
xlabel('Residual Norm || Gm-d ||_2');
ylabel('Solution Norm || m ||_2');

disp('Displaying model norm vs. residual norm for CGLS solution (fig. 5)')

%
% Now, use CGLS to solve the explicitly regularized problem for several
% values of alpha.
%

disp('Note: The following computations may take several minutes on most computers.');
alphas=[1.0e-6; 3.0e-6; 5.0e-6; 7.0e-6;  1.0e-5; 3.0e-5; 5.0e-5; ...
        6.0e-5; 7.0e-5; 8.0e-5; 9.0e-5; 1.0e-4; 2.0e-4;  3.0e-4; 4.0e-4; ...
        5.0e-4; 6.0e-4; 7.0e-4; 1.0e-3; 1.5e-3; 2.0e-3; 2.5e-3];


% Get size of G
[m,n]=size(G);

% Initialize rho, eta and solution matrices for explicit calcluations
explicitrho=zeros(length(alphas),1);
expliciteta=zeros(length(alphas),1);
explicitsols=zeros(40000,length(alphas));

% Use CGLS to get solutions at alpha values listed in vector above
for i=1:length(alphas)
  Gtemp=[G; alphas(i)*speye(n,n)];
  dtemp=[dn; zeros(n,1)];
  [X,rho2,eta2]=cgls(Gtemp,dtemp,200);
  explicitrho(i)=norm(G*X(:,200)-dn);
  expliciteta(i)=norm(X(:,200));
  explicitsols(:,i)=X(:,200);
end

% Plot both L-curves
figure(6)
clf
loglog(rho(1:200),eta(1:200),'ko-');
hold on
loglog(explicitrho,expliciteta,'kx-');
xlabel('Residual Norm || Gm-d ||_2');
ylabel('Solution Norm || m ||_2');
bookfonts
legend('CGLS, no regularization','CGLS, explicit regularization');

disp('Displaying L-curves for CGLS deblurring (fig. 6)')
print -deps c6fblurlcurve.eps 

%
% The corner of the L-curve for the explicitly reguarlized solution is
% at alpha=7e-4, eta=2.77e+4, rho=2.72.
%
% Plot this solution
figure(7)
clf
mexpl=explicitsols(:,18);
imagesc(reshape(mexpl,200,200));
colormap(gray);
H=gca;
set(H,'XTick',[]);
set(H,'YTick',[]);
bookfonts

disp(['Displaying explicitly regularized CGLS solution for alpha ='...
    '7e-4 (fig. 7)'])
print -deps c6fcglsexpl.eps
