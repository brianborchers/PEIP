% Computes forward model derivatives for the EM-38 problem.
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
function J=jac(sigma)

% Number of layers.
M=11;

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
DELTA=sqrt(2/(sigma(1)*MU0*OMEGA));

% Heights of measurements.
heights=[0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.75; 1.0; 1.50];

% Now, do the predictions.
pred=[];
pred=zeros(length(heights),2);
for i=1:length(heights)
  H=heights(i);
  [predv,predh]=predict(R,DELTA,H,M,MU,MU0,sigma,SIGMA0,D,OMEGA);
  pred(i,:)=[predv predh];
end

% Stack up the predictions (vertical over horizontal)
f=[pred(:,1); pred(:,2)];

% Now, loop through the model parameters, computing derivatives.
J=zeros(length(f),M);
I=eye(M);
h=1.0e-8;
for i=1:M
  J(:,i)=(fun(sigma+h*I(:,i))-fun(sigma))/h;
end
