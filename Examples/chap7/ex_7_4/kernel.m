% Example 7.4
% This is the kernel function for the Skaggs and Kabala problem
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% 
% xx       a vector of positions 
% tt       a vector of times
% params   a 2 element array where the first element is the velocity 
%          and the second is the diffusion coefficient
%
% f        the concentration distance xx from the source with tt time till the 
%          sample
function f=kernel(xx,tt,params)

f=xx./(2*sqrt(pi*params(2)*tt.^3)).*...
	  exp(-(xx-params(1).*tt).^2./(4*params(2)*tt));
