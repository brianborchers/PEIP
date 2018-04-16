% Example 5.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear

%
% Data:
%
%   Mass of the Earth (M_e) and Moment of Inertia of the Earth (I_e)
%   (after rescaling)
    d=[5.973; 8.02];

%   Error tolerances on [M_e; I_e]  
    sigma=[0.0005; 0.005];

%   Radius of the Earth (km)
    Re=6370.8;

% throughout we will let g1 be the shellMass as computed by shellMass.m
%                        g2 be the shellInertia as computed by shellInertia.m
% Get the normalizing constraint by evaluating 
% q_i = int(g_i,r,0,1) for i=1,2 using rescaled 
% values of g_i from from shellMass and shellInertia
% Populate q vector using offline evaluation of integral above
q=[1.083221147; 1.757951474];

% determine the values to use a 1000 km
rad = 1000;
%Find radius of interest as a fraction of R_e
ri=rad/Re;

% Get the H matrix for this radius.  
% The following formulas for the integrals were done with Maple.
H=zeros(2,2);
H(1,1)=1.508616069 - 3.520104161*ri + 2.112062496*ri^2;
H(1,2)=3.173750352 - 7.140938293*ri + 4.080536168*ri^2;
H(2,1)=H(1,2);
H(2,2)=7.023621326 - 15.45196692*ri + 8.584426066*ri^2;

% Get the c coeffcients using quadlin to solve the associated 
% Lagrange multiplier problem.
[c,~]=quadlin(H,q',[1.0]);
fprintf('Coefficients for ri=%f (%f km) are %f, %f\n',[ri,ri*Re,c(1),c(2)]);
fprintf('Estimate of density for ri=%f is %f \n',[ri,c'*d]);
fprintf('Standard deviation of density for ri=%f is %f\n\n',...
    [ri,sqrt(c(1)^2*sigma(1)^2+c(2)^2*sigma(2)^2)]);

% compute the averaging kernel.
r=(0.0:0.02:1.0)';
a1000=c(1)*shellMass(r)+c(2)*shellInertia(r);

% now the values at 5000 km
rad = 5000;
%Find radius of interest as a fraction of R_e
ri=rad/Re;

% Get the H matrix.  The following formulas for the integrals were done
% with Maple.
H=zeros(2,2);
H(1,1)=1.508616069 - 3.520104161*ri + 2.112062496*ri^2;
H(1,2)=3.173750352 - 7.140938293*ri + 4.080536168*ri^2;
H(2,1)=H(1,2);
H(2,2)=7.023621326 - 15.45196692*ri + 8.584426066*ri^2;

% Get the c coeffcients using quadlin to solve the associated 
% Lagrange multiplier problem.
[c,lambda]=quadlin(H,q',[1.0]);
fprintf('Coefficients for ri=%f (%f km) are %f, %f\n',[ri,ri*Re,c(1),c(2)]);
fprintf('Estimate of density for ri=%f is %f \n',[ri,c'*d]);
fprintf('Standard deviation of density for ri=%f is %f\n\n',...
    [ri,sqrt(c(1)^2*sigma(1)^2+c(2)^2*sigma(2)^2)]);

% Plot the averaging kernel.
r=(0.0:0.02:1.0)';
a5000=c(1)*shellMass(r)+c(2)*shellInertia(r);

% now plot the averaging kernels
figure(1)
clf

%plot the averaging for 1000 km
plot(r*Re,a1000,'k');
H1=text(350,0.8,'r=1000 km');
set(H1,'FontSize',18);
hold on

%plot the averaging for 5000 km
plot(r*Re,a5000,'k');
H2=text(3000,0.55,'r=5000 km');
set(H2,'FontSize',18);
xlabel('Radius (km)');
ylabel('Averaging Kernel, A(r)');
bookfonts
print -deps2 akernels.eps
