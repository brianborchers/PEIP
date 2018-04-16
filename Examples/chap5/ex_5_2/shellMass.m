% Calculates the value of 4*Pi*rhat^2 for the rescaled version of
% the spherical Earth model (Example 5.1).  
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
function shellMass=shellMass(r)
%
%
% Input Parameters:
%   r - vector of earth radius values (rescaled)
%
% Output Parameters:
%   g1 - vector of values of 4*Pi*r^2*(rescaling constant)
%        
%
% The rescaling constant 1.0344 is obtained as an approximation of  
% 4*(R_e)^3*1000)/10^24, where (R_e)^3*1000)/10^24 is a
% constant obtained by rescaling the problem with:
%   R_e = 6.3708e6, rhat = r/R_e, rhohat = rho/1000, and
%   M_ehat = M_e/10^24
shellMass=1.0344*pi*r.^2;

