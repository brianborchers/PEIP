% Calculates the value of 8/3*Pi*rhat^4 for the
% rescaled version of the spherical Earth model (Example 5.1).  
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
function shellInertia=shellInertia(r)
%
%
% Input Parameters:
%   r - vector of earth radius values (rescaled)
%
% Output Parameters:
%   shellInertia - vector of values of 8/3*Pi*r^4*(rescaling constant)
%
% The constant 2.797866667 is obtained as an approximation of  
% 8/3*(R_e)^5*1000)/10^37, where (R_e)^5*1000)/10^37 is a
% constant obtained by rescaling the problem with:
%   R_e = 6.3708e6, rhat = r/R_e, rhohat = rho/1000, and
%   I_ehat = I_e/10^37
shellInertia=2.797866667*pi*r.^4;

