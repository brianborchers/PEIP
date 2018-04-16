% Example 11.3 slowness function calculation
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
function s = slowfunctionsmooth(depth)
% s = slowfunctionsmooth(depth) 
%
% Returns a vector of true one-dimensional slownesses
%
% Input Parameters:
%   depth - vector of depth values (m)
%
% Output parameters:
%   s - vector of slowness vector (s/km)

% Get velocity vector
v=(3.0+sqrt(depth)/sqrt(1000))*1000;

% Get true 1-d slowness structure
s=1./v;
