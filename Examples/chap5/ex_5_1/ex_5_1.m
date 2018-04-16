% Examples 5.1
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear
rand('seed',0);
randn('seed',0);

%Example 5.1
%We will solve the Shaw problem using a set of basis set of  sines and cosines, plus
%the constant function, on the -pi/2 to p/2 interval at 100 equally spaced
%points for 20 data points

%dicretization of model space
npts=1000;
dth=pi/npts;
th=linspace(-pi/2+dth/2,pi/2-dth/2,npts)';

%Build the basis functions and store them as columns of M
%model point spacing
%Number of sin and cos basis functions to use
%(including the constant function, M(1,:))
T=pi;
%number of basis functions
n=31;
M(:,1)=ones(npts,1);
for i=2:(n+1)/2
    M(:,i)=sin(2*pi*(i-1)*th/T);
end
for i=(n+1)/2+1:n
    M(:,i)=cos(2*pi*(i-(n+1/2))*th/T);
end
%normalize each model basis function to equal length in the appropriate L2
%space
for i=1:n
    M(:,i)=M(:,i)/norm(M(:,i));
end

figure(1)
clf
plot(th,M)
xlabel('\theta (radians)')
ylabel('Basis Function Aplitude')
title('Fourier Basis')
bookfonts

%data point spacing
m=20;
ds=pi/m;
s=linspace(-pi/2+ds/2,pi/2-ds/2,m)';

%construct the G matrix; G_ij is the Shaw forward problem for the ith data
%point and jth model basis function, so each column of G is composed of the
%predicted data for corresponding basis function.

for i=1:n
    Gb(:,i)=shawforward(s,M(:,i),th,dth);
end

%generate set of test models and data
%test model is a triangle function of width nwin
nwin=199;
hwin=(nwin-1)/2;
for i=150:50:850
mtrue=zeros(npts,1);
%we'll use a simple true modelfunction centered on model element
w = 2*(0:(nwin-1)/2)/(nwin-1);
bwin = [w w((nwin-1)/2:-1:1)]';
mtrue(i-hwin:i+hwin)=bwin;
%Comment above, and uncomment below for a spike model test
%mtrue(i)=1.;
dtrue=shawforward(s,mtrue,th,dth);
dn=dtrue+(1e-6)*randn(size(dtrue));

%discrepancy principle residual norm target
disc=sqrt(m)*1e-6;

%search the L-curve 
alphas=10.^(-6:.05:1)';
for j=1:length(alphas)
    alpha=alphas(j);
Gba=[Gb ; alpha*eye(n,n)];
beta=Gba\[dn;zeros(n,1)];
chi=norm(Gb*beta-dn);

%do a check to see if we have crossed the discrepancy principle defined
%misfit value
if chi-disc > 0
    %yes, average with the previous alpha and quit the loop
    if j > 1
    alpha=(alpha+alphas(j-1))/2;
    end
    break
end
end

%calulate a solution for beta and the corresponding 2-norm misfit for this
%particulary alpha value
Gba=[Gb; alpha*eye(n,n)];
beta=Gba\[dn;zeros(n,1)];
chi=norm(Gb*beta-dn);

%build the model
mod=zeros(size(M(:,1)));
for j=1:n
    mod=mod+beta(j)*M(:,j);
end

figure(2)
clf
plotconst(mod,-pi/2,pi/2);
hold on
H=plotconst(mtrue,-pi/2,pi/2);
set(H,'linestyle','--')
legend('Recovered Model','True Model','location','northwest')
hold off
ylabel('Intensity')
xlabel('\theta')
bookfonts

disp(['alpha, chi = ',num2str(alpha),', ',num2str(chi)]);
ylim([-0.55 1])
pause(0.2)
end
print -deps2 Shaw_Fourier.eps


