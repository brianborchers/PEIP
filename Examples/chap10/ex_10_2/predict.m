% EM forward model function
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% [predv,predh]=predict(R,DELTA,H,M,MU,MU0,SIGMA,SIGMA0,D,OMEGA)
%
% Given a physical model, predict the EM reading at a height H above the
% ground.  
%
%  R         Coil spacing (m)
%  DELTA     Scaling factor.  Typically sqrt(2/SIGMA(1)*MU0*OMEGA)
%  H         height above ground (m)
%  M         Number of layers.
%  MU        Array of magnetic permeabilities per layer.
%  MU0       Magnetic permeability of the air above ground.
%  SIGMA     Array of conductivities per layer
%  SIGMA0    Conductivity of the air above ground (0)
%  D         Array of layer thicknesses (m)
%  OMEGA     Angular frequency (2*pi*f) (rad/s)
%
%  predv     Predicted reading with coils oriented vertically
%  predh     Predicted reading with coils oriented horizontally
%
function [predv,predh]=predict(R,DELTA,H,M,MU,MU0,SIGMA,SIGMA0,D,OMEGA)
[WT0,WT1]=hankelwts;
g=hankelpts(R/DELTA);
R0=r0(g/DELTA,R,DELTA,H,M,MU,MU0,SIGMA,SIGMA0,D,OMEGA);
F0=-R0.*g.*g.*exp(-2*g*H/DELTA);
F1=-R0.*g.*exp(-2*g*H/DELTA);
T0=WT0.'*F0/(R/DELTA);
T2=WT1.'*F1/(R/DELTA);
predv=imag(1+T0*(R/DELTA)^3)*1000*4/(MU0*OMEGA*R*R);
predh=imag(1+T2*(R/DELTA)^2)*1000*4/(MU0*OMEGA*R*R);

