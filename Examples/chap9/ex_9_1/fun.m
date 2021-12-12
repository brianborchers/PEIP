% Example 9.1
% Computes the weighted residual for the model p
% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% fvec=fun(p)
%
%
% p is expected to be a two element vector
%   - the first element is transmissivity 
%   - the second element is the storage coefficient
%
function fvec=fun(p)
% global variables, these are 
% H     - the recorded head for each time
% TM    - the times the head was recorded
% SIGMA - the standard deviation for a time
% D     - the distance between the wells
% Q     - the volume of the slug
global H;
global TM;
global SIGMA;
global D;
global Q;

% Compute the function values.
fvec=zeros(length(TM),1);
for i=1:length(TM)
  fvec(i)=(Q*exp(-D^2*p(1)/(4*p(2)*TM(i)))/(4*pi*p(2)*TM(i)) - H(i))/SIGMA(i);
end  
