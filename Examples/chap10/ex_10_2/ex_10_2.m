% Example 10.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% This script creates the synthetic data set for the EM-38 example
%

% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

% Number of layers.
M=11;

% The model conductivities in S/m.
mtrue=[100; 95; 90; 95; 100; 130; 160; 200; 250; 310; 360]/1000.0;

% Layer thicknesses.
D=0.2*ones(10,1);

% Operating frequency.
f=14600;
OMEGA=2*pi*f;

% Distance between coils.
R=1.0;

% Conductivity of the air above the top layer.
SIGMA0=0;

% Magnetic permeabilities.
MU0=pi*4e-7;
MU=MU0*ones(11,1);

% Scaling factor Delta.
DELTA=sqrt(2/(mtrue(1)*MU0*OMEGA));

% Heights of measurements.
heights=[0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.75; 1.0; 1.50];

% Now, do the predictions, using a black box routine
pred=zeros(length(heights),2);
for i=1:length(heights)
  H=heights(i);
  [predv,predh]=predict(R,DELTA,H,M,MU,MU0,mtrue,SIGMA0,D,OMEGA);
  pred(i,:)=[predv predh];
end

% Add random noise.
rand('state',0);
datanf=[pred(:,1); pred(:,2)];
data=datanf+0.1*randn(size(datanf));
%save the data
save EMdata.mat data

% Now, solve the problem
m0=200*ones(11,1)/1000;

% First, try using LM.
global DATA;
DATA=data;

mlm=lm('fund','jac',m0,1.0e-6,50);
disp(['The chi^2 value is ' num2str(norm(fund(mlm),2)^2) ' with 18-11=7 dof']);
Jlm=jac(mlm);
disp(['cond(J''*J) is ' num2str(cond(Jlm'*Jlm),'%1.5e')]);

% Plot the LM solution.
figure(1)
clf

plotconst(mlm*1000,0,2.2);
xlabel('depth (m)');
ylabel('Conductivity (mS/m)');
bookfonts

disp('Displaying the LM solution (fig. 1)')
print -deps c9fmlm.eps


% Next, Occam.
L2=get_l_rough(11,2);
m0=200*ones(11,1)/1000;
delta=0.1*sqrt(18);
moccam=occam('fun','jac',L2,data,m0,delta);

% plot 
figure(2)
clf

plotconst(moccam*1000,0,2.2);
xlabel('depth (m)');
ylabel('Conductivity (mS/m)');
bookfonts

disp('Displaying the Occam''s inversion recovered model (fig. 2)')
print -deps c9fmoccam.eps

% Plot the true model.
figure(3)
clf

plotconst(mtrue*1000,0,2.2);
xlabel('depth (m)');
ylabel('Conductivity (mS/m)');
bookfonts

disp('Displaying the true model (fig. 3)')
print -deps c9fmtrue.eps
